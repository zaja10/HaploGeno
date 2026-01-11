# Script: 06_structure_analysis.R
# Description: Performs Factor Analysis on Local GEBVs to identify Latent Genomic Gradients (Haplo-FA).
#              Can use Covariance (magnitude-sensitive) or Correlation (pattern-sensitive) matrices.
# Dependencies: stats, MASS

library(stats)
# MASS is suggested for ginv(), we call it via namespace.

#' Analyze Block Structure (Haplo-FA)
#'
#' @param local_gebv_mat (N x M) Matrix of Local GEBVs (or data.frame).
#' @param blocks_df Data frame of blocks (must match columns of local_gebv_mat).
#' @param top_n Number of top variable blocks to use (default: 500).
#' @param factors Number of factors to extract (default: 2).
#' @param use_covariance If TRUE, preserves magnitude (Covariance). If FALSE, standardized (Correlation).
#' @return List containing Loadings, Scores, Uniqueness, and metadata.
analyze_block_structure_func <- function(local_gebv_mat, blocks_df, top_n = 500, factors = 2, use_covariance = FALSE) {
    
    # --- 1. Robust Input Preparation ---
    local_gebv_mat <- as.matrix(local_gebv_mat)
    n_total <- ncol(local_gebv_mat)
    
    # Check alignment
    if (nrow(blocks_df) != n_total) {
        # Warning instead of stop, in case user pre-filtered the matrix
        warning("Number of blocks in blocks_df (", nrow(blocks_df), 
                ") does not match columns in GEBV matrix (", n_total, "). Positions may be misaligned.")
    }

    message("Running Haplo-FA on top ", min(top_n, n_total), " blocks...")

    # --- 2. Select High-Variance Blocks ---
    # We filter by variance to focus on "active" regions and reduce dimensionality
    vars <- apply(local_gebv_mat, 2, var, na.rm = TRUE)
    vars[is.na(vars)] <- 0

    efficient_n <- min(top_n, n_total)
    
    # Validation: DoF Check for Factor Analysis
    # Formula: d = 0.5 * ((p-k)^2 - (p+k))
    # If d < 0, the model is under-identified and will crash.
    # Heuristic: p must be roughly > factors + 1
    if (efficient_n <= (factors + 1)) {
        stop("Error: Number of selected blocks (", efficient_n, ") is too small for ", factors, " factors.")
    }

    top_indices <- order(vars, decreasing = TRUE)[1:efficient_n]
    X_sub <- local_gebv_mat[, top_indices, drop = FALSE]

    # --- 3. Factor Analysis Prep ---
    if (use_covariance) {
        message("Using Covariance Matrix (Magnitude Preserved)...")
        # Center but DO NOT scale
        X_anal <- scale(X_sub, center = TRUE, scale = FALSE)
        # Calculate Covariance
        R <- cov(X_sub, use = "pairwise.complete.obs")
    } else {
        message("Using Correlation Matrix (Standardized Pattern)...")
        # Z-score standardization
        X_anal <- scale(X_sub, center = TRUE, scale = TRUE)
        # Calculate Correlation
        R <- cor(X_sub, use = "pairwise.complete.obs")
    }

    # --- 4. Fit Factor Analysis ---
    # stats::factanal works on correlation matrices by default. 
    # Even if we pass a cov matrix, it converts to cor internally. 
    # We must handle the re-scaling manually afterwards.
    
    # We use tryCatch because Maximum Likelihood estimation can fail to converge
    fa_fit <- tryCatch({
        stats::factanal(covmat = R, factors = factors, rotation = "varimax")
    }, error = function(e) {
        stop("Factor Analysis failed to converge. Try reducing 'factors' or 'top_n'. Error: ", e$message)
    })

    # --- 5. Extract & Scale Loadings ---
    L_rot <- fa_fit$loadings[, , drop = FALSE] # Keep as matrix even if 1 factor
    u2 <- fa_fit$uniquenesses
    
    # Scale factors: SDs (if covariance) or 1s (if correlation)
    scale_factors <- sqrt(diag(R))

    # Rescale Loadings: L_cov = L_cor * SD
    L_rot <- sweep(L_rot, 1, scale_factors, "*")
    
    # Rescale Uniqueness: U_cov = U_cor * Var
    u2 <- u2 * (scale_factors^2)

    # --- 6. Compute Factor Scores (Regression Method) ---
    # Formula: F_hat = (X - mu) * R^-1 * L
    # We need to handle NAs in X_anal before multiplication
    
    # Simple Imputation: Replace NA with 0 (Mean) since X_anal is centered
    if (any(is.na(X_anal))) {
        X_anal[is.na(X_anal)] <- 0
    }

    # Calculate weights: R^-1 * L
    # Use Generalized Inverse if R is singular (common in genomic data due to LD)
    weights <- tryCatch({
        solve(R, L_rot)
    }, error = function(e) {
        message("Note: Matrix is singular. Using generalized inverse for scores.")
        MASS::ginv(R) %*% L_rot
    })

    S_rot <- X_anal %*% weights
    colnames(S_rot) <- paste0("Factor", 1:factors)

    # --- 7. Package Results ---
    # Get positions if available
    if ("Start" %in% names(blocks_df) && "End" %in% names(blocks_df)) {
        # Map indices to block rows
        # If blocks_df matches total cols, we use top_indices directly
        blk_pos <- (blocks_df$Start[top_indices] + blocks_df$End[top_indices]) / 2
    } else {
        blk_pos <- top_indices
    }

    return(list(
        Loadings = L_rot,
        Scores = S_rot,
        Uniqueness = u2,
        BlockIndices = top_indices,
        Pos = blk_pos,
        VarExplained = colSums(L_rot^2) # Approximate variance explained per factor
    ))
}
