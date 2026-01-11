# Script: 07_vis_genotype_structure.R
# Description: Visualization functions for Genotype and Haplotype Structure (Blocks, HRM, PCA).
# Dependencies: graphics (base R)

#' Plot Distribution of Haplotype Block Sizes
#'
#' @param blocks_df Data frame with columns Start and End.
#' @param main Title of the plot.
#' @param col Color of histogram bars.
plot_block_sizes <- function(blocks_df, main = "Distribution of Haplotype Block Sizes", col = "steelblue") {
    if (missing(blocks_df) || is.null(blocks_df)) stop("blocks_df is required.")

    sizes <- blocks_df$End - blocks_df$Start + 1

    hist(sizes,
        breaks = 30,
        col = col,
        border = "white",
        main = main,
        xlab = "Number of Markers per Block"
    )

    # Add median line
    med_val <- median(sizes)
    abline(v = med_val, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median:", med_val), col = "red", lty = 2, bty = "n")
}

#' Plot Haplotype Relationship Matrix (HRM)
#'
#' @param hrm_mat Square numeric matrix (HRM).
#' @param main Title of the plot.
plot_hrm <- function(hrm_mat, main = "Haplotype Relationship Matrix") {
    if (missing(hrm_mat) || is.null(hrm_mat)) stop("hrm_mat is required.")

    # Heatmap palette (Yellow-Red for relatedness)
    pal <- colorRampPalette(c("white", "orange", "darkred"))(100)

    # Image plot (flip Y to put sample 1 at top left visually)
    image(1:nrow(hrm_mat), 1:ncol(hrm_mat), t(hrm_mat)[, nrow(hrm_mat):1],
        col = pal,
        axes = FALSE,
        main = main,
        xlab = "Individuals",
        ylab = "Individuals"
    )

    box()
}

#' Plot PCA of Haplotype Relationship Matrix
#'
#' @param hrm_mat Square numeric matrix (HRM).
#' @param groups Optional vector of groups for coloring (length = nrow(hrm_mat)).
#' @param main Title of the plot.
plot_pca <- function(hrm_mat, groups = NULL, main = "PCA of Haplotype Relationship Matrix") {
    if (missing(hrm_mat) || is.null(hrm_mat)) stop("hrm_mat is required.")

    message("Calculating Eigen decomposition of HRM...")
    eig <- eigen(hrm_mat, symmetric = TRUE)

    # Extract top 2 PCs
    pc1 <- eig$vectors[, 1]
    pc2 <- eig$vectors[, 2]

    # Variance explained
    var_expl <- eig$values / sum(eig$values) * 100

    # Plot
    xlab <- paste0("PC1 (", round(var_expl[1], 1), "%)")
    ylab <- paste0("PC2 (", round(var_expl[2], 1), "%)")

    col_vec <- "black"
    if (!is.null(groups)) {
        if (length(groups) != nrow(hrm_mat)) {
            warning("Length of groups does not match number of individuals. Ignoring groups.")
        } else {
            # Convert groups to factor to get integer codes for colors
            groups <- as.factor(groups)
            col_vec <- as.integer(groups)
        }
    }

    plot(pc1, pc2,
        xlab = xlab, ylab = ylab,
        main = main,
        pch = 19,
        col = col_vec
    )

    if (!is.null(groups)) {
        legend("topright", legend = levels(groups), col = seq_along(levels(groups)), pch = 19)
    }
}
