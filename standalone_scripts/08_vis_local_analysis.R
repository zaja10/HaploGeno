# Script: 08_vis_local_analysis.R
# Description: Visualization functions for Local Analyses (GEBVs, Factor Analysis).
# Dependencies: graphics (base R)

#' Plot Local GEBV Heatmap
#'
#' @param local_gebv_mat Matrix of Local GEBVs (Individuals x Blocks).
#' @param block_range Optional integer vector specifying which blocks (columns) to plot.
#' @param main Title of the plot.
plot_gebv_image <- function(local_gebv_mat, block_range = NULL, main = "Local GEBV Heatmap") {
    if (missing(local_gebv_mat) || is.null(local_gebv_mat)) stop("local_gebv_mat is required.")

    mat <- local_gebv_mat

    if (!is.null(block_range)) {
        if (max(block_range) > ncol(mat)) stop("Block range exceeds number of blocks.")
        mat <- mat[, block_range, drop = FALSE]
    }

    message("Plotting GEBV heatmap...")

    # Palette
    pal <- colorRampPalette(c("blue", "white", "red"))(100)

    # Transpose for image() which expects x as rows, y as cols
    # We want X-axis = Blocks, Y-axis = Individuals

    image(
        x = 1:ncol(mat),
        y = 1:nrow(mat),
        z = t(mat),
        col = pal,
        xlab = "Block Index",
        ylab = "Individual Index",
        main = main,
        useRaster = TRUE
    )
}

#' Plot Factor Analysis Heatmap (Genetic Correlation)
#'
#' Visualization of the model-implied correlation matrix from Haplo-FA.
#'
#' @param fa_results List object returned by `analyze_block_structure_func`.
#' @param subset_indices Optional vector of block indices to zoom in on.
#' @param show_values Logical. Overlay correlation values?
#' @param label_axes Logical. Label axes with Block IDs?
#' @param main Title of the plot.
#' @param ... Additional arguments passed to image().
plot_factor_heatmap <- function(fa_results, subset_indices = NULL, show_values = FALSE, label_axes = TRUE, main = "Latent Block Correlations", ...) {
    if (missing(fa_results) || is.null(fa_results)) stop("fa_results is required.")

    L <- fa_results$Loadings
    Psi <- fa_results$Uniqueness # Note: Standardized uniqueness scaled by Var
    block_ids_all <- fa_results$BlockIndices

    # 1. Reconstruct Correlation Matrix
    # Sigma = L L' + Psi (Diagonal)
    Sigma <- tcrossprod(L)
    diag(Sigma) <- diag(Sigma) + Psi

    # Convert to Correlation
    G_cor <- cov2cor(Sigma)

    # 2. Handle Subsetting
    if (!is.null(subset_indices)) {
        match_idx <- match(subset_indices, block_ids_all)
        valid <- !is.na(match_idx)

        if (sum(valid) < 2) stop("Fewer than 2 specified blocks were found in the FA results.")

        target_idx <- match_idx[valid]
        G_plot <- G_cor[target_idx, target_idx, drop = FALSE]
        labels <- block_ids_all[target_idx]
    } else {
        G_plot <- G_cor
        labels <- block_ids_all
    }

    # 3. Plotting
    n <- nrow(G_plot)
    pal <- colorRampPalette(c("navy", "white", "firebrick"))(100)

    image(1:n, 1:n, t(G_plot)[, n:1],
        col = pal, axes = FALSE,
        xlab = "Block Index", ylab = "Block Index",
        main = main, ...
    )

    # 4. Show Values
    if (show_values) {
        grid_x <- rep(1:n, each = n)
        grid_y <- rep(n:1, times = n)
        val_vec <- as.vector(t(G_plot))

        text_col <- ifelse(abs(val_vec) > 0.5, "white", "black")

        text(
            x = grid_x, y = grid_y, labels = round(val_vec, 2),
            col = text_col, cex = 0.7, font = 1
        )
    }

    # 5. Axes
    if (label_axes) {
        if (!is.null(subset_indices) || n <= 20) {
            at_idx <- 1:n
        } else {
            at_idx <- seq(1, n, floor(n / 20))
            labels <- labels[at_idx]
        }

        axis(1, at = at_idx, labels = labels, las = 2, cex.axis = 0.8)

        # Correct Y-axis labels (flipped orientation)
        y_at <- n - at_idx + 1
        axis(2, at = y_at, labels = labels, las = 2, cex.axis = 0.8)
    }

    box()
}
