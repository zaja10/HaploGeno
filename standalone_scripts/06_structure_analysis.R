# Script: 06_structure_analysis.R
# Description: Factor Analysis for Latent Genomic Gradients (Haplo-FA).
# Dependencies: stats, MASS

#' Analyze Block Structure (Haplo-FA)
#'
#' @param local_gebv_mat (N x M) Matrix of Local GEBVs.
#' @param blocks_df Data frame of blocks.
#' @param top_n Number of top variable blocks to use.
#' @param factors Number of factors.
#' @param use_covariance If TRUE, uses cov() instead of cor().
#' @return List containing Loadings, Scores, Communality.
analyze_block_structure_func <- function(local_gebv_mat, blocks_df, top_n = 500, factors = 2, use_covariance = FALSE) {
    message("Running Haplo-FA on top ", top_n, " blocks...")

    # 1. Select High-Variance Blocks
    vars <- apply(local_gebv_mat, 2, var)
    vars[is.na(vars)] <- 0

    n_total <- ncol(local_gebv_mat)
    efficient_n <- min(top_n, n_total)

    top_indices <- order(vars, decreasing = TRUE)[1:efficient_n]
    X_sub <- local_gebv_mat[, top_indices, drop = FALSE]

    # 2. Factor Analysis Prep
    if (use_covariance) {
        # Centered but not scaled
        X_anal <- scale(X_sub, center = TRUE, scale = FALSE)
        R <- cov(X_sub, use = "pairwise.complete.obs")
    } else {
        # Z-scores
        X_anal <- scale(X_sub, center = TRUE, scale = TRUE)
        R <- cor(X_sub, use = "pairwise.complete.obs")
    }

    # 3. Fit Factor Analysis
    # using stats::factanal with rotation
    fa_fit <- stats::factanal(covmat = R, factors = factors, rotation = "varimax")

    # 4. Extract & Scale Loadings
    L_rot <- fa_fit$loadings[, ]
    u2 <- fa_fit$uniquenesses
    scale_factors <- sqrt(diag(R))

    # Rescale back to original variance units
    L_rot <- sweep(L_rot, 1, scale_factors, "*")
    u2 <- u2 * (scale_factors^2)

    # 5. Compute Factor Scores (Regresison Method)
    # Scores = X_anal * R^-1 * L
    weights <- tryCatch(solve(R, L_rot), error = function(e) MASS::ginv(R) %*% L_rot)
    S_rot <- X_anal %*% weights

    return(list(
        Loadings = L_rot,
        Scores = S_rot,
        Uniqueness = u2,
        BlockIndices = top_indices,
        Pos = (blocks_df$Start[top_indices] + blocks_df$End[top_indices]) / 2
    ))
}
