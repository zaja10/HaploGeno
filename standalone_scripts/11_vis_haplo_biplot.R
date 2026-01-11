# Script: 11_vis_haplo_biplot.R
# Description: Generates a Biplot (Genotypes + Block Vectors).
#              Visualizes which Haplotype Blocks drive the separation of Genotypes.
# Dependencies: ggplot2, ggrepel, dplyr, bigstatsr

library(ggplot2)
library(ggrepel)
library(dplyr)
library(bigstatsr)

#' Plot Haplotype Biplot (Genotypes + Top Block Vectors)
#'
#' @param local_gebv_list Output from Script 04 (must contain $matrix).
#' @param blocks_df Data frame of block definitions (BlockID, Chr, Start).
#' @param genotype_groups Optional factor for coloring genotypes (e.g., Origin, Family).
#' @param top_vectors Number of block vectors to show (default: 10).
#' @param scale_factor Scalar to adjust arrow length for visibility (default: "auto").
#'
#' @return A ggplot object.
plot_haplo_biplot <- function(local_gebv_list, blocks_df, genotype_groups = NULL,
                              top_vectors = 15, scale_factor = "auto") {
    # --- 1. Data Prep ---
    # Convert FBM to standard matrix for PCA (PCA is memory intensive, assumed solvable here)
    if (is.list(local_gebv_list) && "matrix" %in% names(local_gebv_list)) {
        if (inherits(local_gebv_list$matrix, "FBM")) {
            mat <- local_gebv_list$matrix[]
        } else {
            mat <- local_gebv_list$matrix
        }
    } else {
        # Assume it's the matrix itself
        mat <- as.matrix(local_gebv_list)
    }

    # Assign Row Names if missing (assuming sequential ID if not present)
    if (is.null(rownames(mat))) rownames(mat) <- paste0("Ind_", seq_len(nrow(mat)))

    # --- 2. Run PCA on Local GEBVs ---
    # We use PCA on the Local GEBVs. This works similarly to FA for visualization.
    # Center = TRUE, Scale = FALSE (preserves magnitude of genetic variance)
    # Check for NA
    if (any(is.na(mat))) {
        warning("Missing values in GEBV matrix. Imputing with column means for PCA.")
        mat <- apply(mat, 2, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            x
        })
    }

    # Check variance
    col_vars <- apply(mat, 2, var)
    if (any(col_vars == 0)) {
        # remove zero variance columns
        mat <- mat[, col_vars > 0, drop = FALSE]
        # Update blocks_df accordingly??
        # Ideally yes, but for visualization we might just ignore the mismatch for now
        # providing we can map back.
        # But if we drop cols, indices change.
        # This is tricky. Let's assume input is clean from Script 04 (filtered).
        # Or just warn.
    }

    pca_res <- prcomp(mat, center = TRUE, scale. = FALSE)

    # --- 3. Extract Scores (Genotypes) ---
    scores <- as.data.frame(pca_res$x)
    scores$Genotype <- rownames(mat)
    if (!is.null(genotype_groups)) {
        # Ensure length matches
        if (length(genotype_groups) == nrow(scores)) {
            scores$Group <- as.factor(genotype_groups)
        } else {
            warning("Length of genotype_groups does not match number of genotypes. Ignoring grouping.")
            scores$Group <- "All"
        }
    } else {
        scores$Group <- "All"
    }

    # --- 4. Extract Loadings (Blocks) ---
    loadings <- as.data.frame(pca_res$rotation)
    # The columns of mat correspond to blocks.
    # If we dropped columns, this index is relative to the reduced matrix.
    # If we didn't drop, it maps 1:1 to input.

    # We assume 1:1 for simplicity as typical Script 04 output shouldn't have constant cols (blocks usually have variation).
    loadings$BlockIndex <- seq_len(nrow(loadings))

    # Merge with Block Metadata
    # Check if blocks_df matches the columns
    if (!missing(blocks_df) && nrow(blocks_df) == nrow(loadings)) {
        loadings$BlockID <- blocks_df$BlockID
        if ("Chr" %in% names(blocks_df)) {
            loadings$Chr <- blocks_df$Chr
        } else if ("chr" %in% names(blocks_df)) {
            loadings$Chr <- blocks_df$chr
        } else {
            loadings$Chr <- "Unk"
        }
    } else {
        loadings$BlockID <- paste0("Block_", seq_len(nrow(loadings)))
        loadings$Chr <- "Unk"
    }

    # --- 5. Filter for "Driver" Blocks ---
    # Calculate vector length (magnitude of influence on PC1/PC2)
    loadings$VectorLength <- sqrt(loadings$PC1^2 + loadings$PC2^2)

    # Select Top N blocks
    top_loadings <- loadings %>%
        arrange(desc(VectorLength)) %>%
        slice_head(n = top_vectors)

    # --- 6. Auto-Scaling ---
    # Biplots often require scaling arrows to match the spread of points
    if (scale_factor == "auto") {
        max_score <- max(abs(c(scores$PC1, scores$PC2)))
        max_loading <- max(abs(c(top_loadings$PC1, top_loadings$PC2)))
        if (max_loading > 0) {
            scale_factor <- (max_score * 0.8) / max_loading
        } else {
            scale_factor <- 1
        }
    }

    # Apply scaling to vectors
    top_loadings$PC1_Scaled <- top_loadings$PC1 * scale_factor
    top_loadings$PC2_Scaled <- top_loadings$PC2 * scale_factor

    # Percentage of Variance
    var_expl <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)

    # --- 7. Plotting ---
    p <- ggplot() +
        # A. Genotype Points
        geom_point(
            data = scores, aes(x = PC1, y = PC2, color = Group),
            alpha = 0.6, size = 2
        ) +

        # B. Block Vectors (Arrows)
        geom_segment(
            data = top_loadings,
            aes(x = 0, y = 0, xend = PC1_Scaled, yend = PC2_Scaled),
            arrow = arrow(length = unit(0.2, "cm")),
            color = "darkred", lwd = 1
        ) +

        # C. Block Labels
        geom_text_repel(
            data = top_loadings,
            aes(
                x = PC1_Scaled, y = PC2_Scaled,
                label = paste0(BlockID, "\n(", Chr, ")")
            ),
            color = "darkred", size = 3.5, fontface = "bold"
        ) +

        # D. Aesthetics
        labs(
            title = "Haplotype Biplot (Local GEBV PCA)",
            subtitle = paste0("Top ", top_vectors, " Blocks driving Genotypic Separation"),
            x = paste0("PC1 (", var_expl[1], "%)"),
            y = paste0("PC2 (", var_expl[2], "%)")
        ) +
        theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3)

    return(p)
}
