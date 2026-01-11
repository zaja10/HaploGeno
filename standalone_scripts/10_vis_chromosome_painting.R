# Script: 10_vis_chromosome_painting.R
# Description: Visualizes Local GEBV distribution along chromosomes (Chromosome Painting).
# Dependencies: ggplot2, dplyr, ggrepel

library(ggplot2)
library(dplyr)
library(ggrepel)

#' Plot Chromosome Painting of Local GEBVs
#'
#' Visualizes the composition of genetic merit along a specific chromosome for selected genotypes.
#'
#' @param local_gebv_list List containing $matrix (GEBVs) or a matrix object itself.
#' @param blocks_df Data frame with BlockID, Chr, Start, End.
#' @param total_gebv_vec (Optional) Vector of total GEBVs for sorting/labeling.
#' @param significance_df (Optional) Data frame with BlockID and P_Value for highlighting.
#' @param target_chr The chromosome to plot (e.g., "1A").
#' @param target_genotypes Vector of genotype IDs to plot.
#' @param highlight_genes_df (Optional) Data frame with Gene, Chr, Pos, Trait_Category.
#' @param map_data (Optional) Map data to map indices to positions. If NULL, uses indices.
#' @return A ggplot object.
#' @export
plot_chromosome_painting <- function(local_gebv_list, blocks_df, total_gebv_vec = NULL,
                                     significance_df = NULL, target_chr, target_genotypes,
                                     highlight_genes_df = NULL, map_data = NULL) {
    # 1. Unpack GEBVs
    if (is.list(local_gebv_list) && "matrix" %in% names(local_gebv_list)) {
        gebv_mat <- local_gebv_list$matrix
    } else {
        gebv_mat <- as.matrix(local_gebv_list)
    }

    # 2. Filter Blocks for Target Chromosome
    if (!"Chr" %in% names(blocks_df)) stop("blocks_df must contain 'Chr' column.")

    # Ensure strict string matching for Chr
    chr_blocks <- blocks_df[as.character(blocks_df$Chr) == as.character(target_chr), ]

    if (nrow(chr_blocks) == 0) {
        warning(paste("No blocks found for Chromosome", target_chr))
        return(NULL)
    }

    # 3. Handle Physical Positions vs Indices
    # If map_data is provided, map Start/End indices to Positions
    # Otherwise check if blocks_df already has positions?
    # Usually blocks_df from Script 01 has indices.

    plot_data <- chr_blocks
    xlabel <- "Marker Index"

    if (!is.null(map_data)) {
        # Check if map_data has positions
        pos_col <- intersect(names(map_data), c("pos", "bp", "position", "start"))
        if (length(pos_col) > 0) {
            # Map indices to positions
            pvec <- map_data[[pos_col[1]]]
            plot_data$StartPos <- pvec[plot_data$Start]
            plot_data$EndPos <- pvec[plot_data$End]
            xlabel <- paste("Position (", pos_col[1], ")", sep = "")
        } else {
            # Fallback to indices
            plot_data$StartPos <- plot_data$Start
            plot_data$EndPos <- plot_data$End
        }
    } else {
        # Fallback to indices
        plot_data$StartPos <- plot_data$Start
        plot_data$EndPos <- plot_data$End
    }

    # 4. Prepare Data for selected genotypes
    res_list <- list()

    available_genos <- intersect(target_genotypes, rownames(gebv_mat))
    if (length(available_genos) == 0) stop("No target_genotypes found in GEBV matrix.")

    for (gid in available_genos) {
        # Extract values for blocks on this chr
        # Assuming gebv_mat columns match BlockIDs (1..M)
        vals <- gebv_mat[gid, chr_blocks$BlockID]

        tmp <- plot_data
        tmp$Genotype <- gid
        tmp$Value <- vals

        # Add Significance info if available
        if (!is.null(significance_df)) {
            # Match by BlockID
            # assuming significance_df has BlockID
            merged <- merge(tmp, significance_df, by = "BlockID", all.x = TRUE)
            # If P_Value exists, flag it
            if ("P_Value" %in% names(merged)) {
                # Top 1% threshold or provided? Let's just store the P-value or a generic flag
                # For painting, maybe just outline significant ones?
                merged$IsSig <- !is.na(merged$P_Value) & merged$P_Value < (0.05 / nrow(blocks_df)) # Bonferroni rough
                res_list[[gid]] <- merged
            } else {
                tmp$IsSig <- FALSE
                res_list[[gid]] <- tmp
            }
        } else {
            tmp$IsSig <- FALSE
            res_list[[gid]] <- tmp
        }
    }

    final_df <- do.call(rbind, res_list)

    # 5. Plotting

    p <- ggplot(final_df, aes(xmin = StartPos, xmax = EndPos, ymin = 0, ymax = 1, fill = Value)) +
        geom_rect() +
        facet_grid(Genotype ~ .) +
        scale_fill_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = 0,
            name = "Local GEBV"
        ) +
        labs(
            title = paste("Chromosome Painting:", target_chr),
            x = xlabel, y = ""
        ) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
        )

    # Add significant borders?
    if (any(final_df$IsSig)) {
        p <- p + geom_rect(
            data = subset(final_df, IsSig),
            aes(xmin = StartPos, xmax = EndPos, ymin = 0, ymax = 1),
            fill = NA, color = "black", linewidth = 0.5, linetype = "dotted"
        )
    }

    # 6. Gene Annotation Overlay
    if (!is.null(highlight_genes_df)) {
        # Filter genes for this Chr
        # Ensure string match
        genes_chr <- highlight_genes_df[as.character(highlight_genes_df$Chr) == as.character(target_chr), ]

        if (nrow(genes_chr) > 0) {
            # genes_df needs Genotype column? No, genes are global per facet
            # BUT ggrepel might struggle with faceting if data doesn't have the facet var.
            # Alternatively, plot genes on a separate track or just overlay on all?
            # Let's overlay on all.
            # To show on all facets, we leave Genotype out of the data passed to geom_text_repel?
            # Yes, if data doesn't have the facet variable, it repeats in all panels.

            p <- p + geom_vline(
                data = genes_chr, aes(xintercept = Pos),
                linetype = "dashed", color = "darkgreen", alpha = 0.5
            ) +
                geom_text_repel(
                    data = genes_chr, aes(x = Pos, y = 1.2, label = Gene),
                    inherit.aes = FALSE, color = "darkgreen", size = 3, force = 5
                )
        }
    }

    return(p)
}
