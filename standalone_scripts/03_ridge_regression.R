# Script: 03_ridge_regression.R
# Description: Estimates marker effects using Ridge Regression (RR-BLUP equivalent).
# Dependencies: bigstatsr, future.apply

library(bigstatsr)
library(future.apply)

#' Estimate Marker Effects
#'
#' @param geno_mat Genotype Matrix or FBM (n x p).
#' @param pheno_vec Phenotype vector (length n).
#' @param lambda Regularization parameter. If "auto", chooses via CV.
#' @param k_folds Number of folds for CV (default: 5).
#' @param n_cores Number of cores for CV.
#' @return A list containing:
#'   - marker_effects: Vector of estimated effects (p x 1).
#'   - model_info: List of model parameters (lambda, intercept, etc).
estimate_marker_effects_func <- function(geno_mat, pheno_vec, lambda = "auto", k_folds = 5, n_cores = 1) {
    n <- nrow(geno_mat)
    p <- ncol(geno_mat)
    y <- as.vector(pheno_vec)

    if (length(y) != n) stop("Phenotype length matches sample size.")

    message("Preparing for Ridge Regression...")

    # 1. Clean Data (Check for Monomorphic Markers)
    # Efficiently check stats
    stats <- bigstatsr::big_colstats(geno_mat)
    keep_bool <- stats$var > 1e-8

    if (!all(keep_bool)) {
        warning("Excluding ", sum(!keep_bool), " monomorphic markers.")
        # We need to act on a subset.
        # For FBM simplicity, we pass ind.col to big_tcrossprodSelf
        ind_col <- which(keep_bool)
    } else {
        ind_col <- bigstatsr::cols_along(geno_mat)
    }

    # 2. Compute GRM (K)
    # K = ZZ' / P  (where Z is centered/scaled genotypes)
    message("Computing GRM...")
    K <- bigstatsr::big_tcrossprodSelf(geno_mat, fun.scaling = bigstatsr::big_scale(), ind.col = ind_col)
    K <- K[] / length(ind_col)

    # 3. Optimize Lambda (if needed)
    if (is.character(lambda) && lambda == "auto") {
        message("Optimizing Lambda via ", k_folds, "-Fold CV...")

        cand_lambdas <- 10^seq(-3, 3, by = 0.5)
        folds <- cut(seq(1, n), breaks = k_folds, labels = FALSE)

        if (n_cores > 1) future::plan(future::multisession, workers = n_cores)

        cv_errors <- future.apply::future_lapply(cand_lambdas, function(lam) {
            err_sum <- 0
            for (k in 1:k_folds) {
                idx_val <- which(folds == k)
                idx_trn <- which(folds != k)

                K_trn <- K[idx_trn, idx_trn]
                K_val_trn <- K[idx_val, idx_trn]
                y_trn <- y[idx_trn]
                y_val <- y[idx_val]

                # Ridge Solution: alpha = (K + lambda*I)^-1 y
                # We assume independent residual variance (I) scaled by lambda relative to K
                # Correct formulation for GBLUP: (G + lambda*I)^-1 y
                alpha <- solve(K_trn + lam * diag(length(idx_trn)), y_trn)

                y_pred <- K_val_trn %*% alpha
                err_sum <- err_sum + sum((y_val - y_pred)^2)
            }
            return(err_sum / n)
        }, future.seed = TRUE)

        best_idx <- which.min(unlist(cv_errors))
        lambda <- cand_lambdas[best_idx]
        message("Optimal Lambda: ", lambda)
    }

    # 4. Final Fit
    message("Fitting final model...")
    I <- diag(n)
    alpha <- solve(K + lambda * I, y)

    # 5. Back-solve for SNP Effects
    # u = Z' alpha
    # u_j = sum_i (Z_ij * alpha_i)
    # Z need to be scaled: (X - center) / scale

    message("Back-calculating marker effects...")

    # Get scaling params
    scaling <- bigstatsr::big_scale()(geno_mat, ind.col = ind_col)
    centers <- scaling$center
    scales <- scaling$scale

    # Raw product X' alpha ("u" unscaled)
    # big_cprodVec calculates: t(X) %*% y
    raw_prod <- bigstatsr::big_cprodVec(geno_mat, alpha, ind.col = ind_col)

    # Adjust for centering/scaling
    # u_hat = (raw_prod - center * sum(alpha)) / scale
    sum_alpha <- sum(alpha)
    u_hat <- (raw_prod - centers * sum_alpha) / scales

    # Expand to full vector including dropped markers (set to 0)
    full_effects <- numeric(p)
    full_effects[ind_col] <- u_hat

    results <- list(
        marker_effects = full_effects,
        model_info = list(
            lambda = lambda,
            intercept = 0, # Should be handled by centering usually
            trait_mean = mean(y)
        )
    )
    return(results)
}
