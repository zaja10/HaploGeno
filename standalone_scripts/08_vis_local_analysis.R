# Script: 08_vis_local_analysis.R
# Description: Visualization functions for Local Analyses (GEBVs, Factor Analysis).
#              Visualizes the genetic architecture and latent structure.
# Dependencies: graphics (base R)

#' Plot Local GEBV Heatmap
#'
#' @param local_gebv_data Matrix of Local GEBVs OR the list returned by calculate_local_gebv_func.
#' @param block_range Optional integer vector specifying which blocks (columns) to plot.
#' @param main Title of the plot.
plot_gebv_image <- function(local_gebv_data, block_range = NULL, main = "Local GEBV Heatmap") {
    if (missing(local_gebv_data) || is.null(local_gebv_data)) stop("local_gebv_data is required.")

    # --- 1. Automatic Unpacking ---
    # Handle if user passed the full list from Script 4
    if (is.list(local_gebv_data) && "matrix" %in% names(local_gebv_data)) {
        mat <- local_gebv_data$matrix
    } else {
        mat <- as.matrix(local_gebv_data)
    }

    # --- 2. Subsetting ---
    if (!is.null(block_range)) {
        if (max(block_range) > ncol(mat)) stop("Error: Block range exceeds number of blocks in matrix.")
        mat <- mat[, block_range, drop = FALSE]
    }

    message("Plotting GEBV heatmap for ", nrow(mat), " individuals x ", ncol(mat), " blocks...")

    # Palette: Blue (Negative) -> White (Zero) -> Red (Positive)
    pal <- colorRampPalette(c("blue", "white", "red"))(100)

    # --- 3. Plotting ---
    # image() expects x=rows, y=cols. To visualize standard matrix orientation:
    # X-axis = Columns (Blocks)
    # Y-axis = Rows (Individuals) - Flipped so Ind 1 is at top
    
    image(
        x = 1:ncol(mat),
        y = 1:nrow(mat),
        z = t(mat)[, nrow(mat):1], # Transpose and flip Y
        col = pal,
        xlab = "Block Index",
        ylab = "Individual Index",
        main = main,
        useRaster = TRUE, # Optimization for large matrices
        axes = FALSE
    )

    # Simple Axes
    axis(1, at = pretty(1:ncol(mat)), cex.axis = 0.8)
    axis(2, at = pretty(1:nrow(mat)), labels = rev(pretty(1:nrow(mat))), cex.axis = 0.8)
    box()
}

#' Plot Factor Analysis Heatmap (Genetic Correlation)
#'
#' Visualization of the model-implied correlation matrix from Haplo-FA.
#' G_smooth = Loadings %*% t(Loadings) + Uniqueness
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
    Psi <- fa_results$Uniqueness 
    block_ids_all <- fa_results$BlockIndices

    # --- 1. Reconstruct Correlation Matrix ---
    # Sigma = L L' + Psi (Diagonal)
    # This represents the "Smoothed" LD structure captured by the factors
    Sigma <- tcrossprod(L)
    
    # Add uniqueness to diagonal (handling vector recycling safely)
    if (length(Psi) != nrow(Sigma)) stop("Dimension mismatch between Loadings and Uniqueness.")
    diag(Sigma) <- diag(Sigma) + Psi

    # Convert to Correlation (Standardizes the covariance matrix)
    G_cor <- cov2cor(Sigma)

    # --- 2. Handle Subsetting ---
    if (!is.null(subset_indices)) {
        # Find positions of requested blocks in the analyzed set
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

    # --- 3. Plotting ---
    n <- nrow(G_plot)
    # Palette: Navy (Neg) -> White (0) -> Firebrick (Pos)
    pal <- colorRampPalette(c("navy", "white", "firebrick"))(100)

    image(1:n, 1:n, t(G_plot)[, n:1],
        col = pal, axes = FALSE,
        xlab = "Block Index", ylab = "Block Index",
        main = main, ...
    )

    # --- 4. Show Values (Optional) ---
    if (show_values) {
        grid_x <- rep(1:n, each = n)
        grid_y <- rep(n:1, times = n) # Y is flipped in image()
        val_vec <- as.vector(t(G_plot))

        # Dynamic text color for contrast
        text_col <- ifelse(abs(val_vec) > 0.5, "white", "black")

        text(
            x = grid_x, y = grid_y, labels = round(val_vec, 2),
            col = text_col, cex = 0.7, font = 1
        )
    }

    # --- 5. Custom Axes ---
    if (label_axes) {
        # Sparse labeling if N is large
        if (!is.null(subset_indices) || n <= 20) {
            at_idx <- 1:n
            lbls <- labels
        } else {
            at_idx <- seq(1, n, length.out = 20)
            lbls <- labels[round(at_idx)]
        }

        # X Axis
        axis(1, at = at_idx, labels = lbls, las = 2, cex.axis = 0.8)

        # Y Axis (Corrected for flip)
        # Data row i is plotted at y-coordinate (n - i + 1)
        y_at <- n - at_idx + 1
        axis(2, at = y_at, labels = lbls, las = 2, cex.axis = 0.8)
    }

    box()
}
