# Script: 03_ridge_regression.R
# Description: Estimates marker effects using Ridge Regression (RR-BLUP).
#              Solves the dual problem: alpha = (K + lambda*I)^-1 y.
# Dependencies: bigstatsr, future, future.apply

library(bigstatsr)
library(future)
library(future.apply)

#' Estimate Marker Effects
#'
#' @param geno_mat Genotype Matrix or FBM (n x p).
#' @param pheno_vec Phenotype vector (length n).
#' @param lambda Regularization parameter. If "auto", chooses via CV.
#' @param k_folds Number of folds for CV (default: 5).
#' @param n_cores Number of cores for CV.
#' @return A list containing:
#'    - marker_effects: Vector of estimated effects (p x 1).
#'    - model_info: List of model parameters (lambda, trait_mean, etc).
estimate_marker_effects_func <- function(geno_mat, pheno_vec, lambda = "auto", k_folds = 5, n_cores = 1) {
    
    # --- 1. Robust Input Handling ---
    # bigstatsr functions require an FBM. If a standard matrix is passed, convert it.
    if (!inherits(geno_mat, "FBM")) {
        message("Converting genotype matrix to FBM format for efficient computation...")
        geno_mat <- bigstatsr::as_FBM(geno_mat)
    }

    n <- nrow(geno_mat)
    p <- ncol(geno_mat)
    y <- as.vector(pheno_vec)

    # Sanity Checks
    if (length(y) != n) stop("Error: Phenotype length does not match sample size (rows of geno_mat).")
    if (any(is.na(y))) stop("Error: Phenotypes contain missing values (NAs). Please impute or filter before analysis.")

    # --- 2. Phenotype Centering ---
    # Ridge Regression shrinks coefficients towards zero. 
    # If y is not centered, the model tries to shrink the intercept, leading to bias.
    y_mean <- mean(y)
    y_centered <- y - y_mean

    message("Preparing for Ridge Regression (N=", n, ", P=", p, ")...")

    # --- 3. Filter Monomorphic Markers ---
    # Efficiently check statistics using bigstatsr
    stats <- bigstatsr::big_colstats(geno_mat)
    keep_bool <- stats$var > 1e-8
    
    if (!all(keep_bool)) {
        n_drop <- sum(!keep_bool)
        warning("Excluding ", n_drop, " monomorphic markers from analysis.")
        ind_col <- which(keep_bool)
    } else {
        ind_col <- bigstatsr::cols_along(geno_mat)
    }

    # --- 4. Compute GRM (K) ---
    # K = ZZ' / P  (where Z is centered/scaled genotypes)
    message("Computing Genomic Relationship Matrix (GRM)...")
    
    # big_tcrossprodSelf computes X %*% t(X) efficiently
    # fun.scaling handles the (X - mean)/sd transformation implicitly
    K <- bigstatsr::big_tcrossprodSelf(geno_mat, fun.scaling = bigstatsr::big_scale(), ind.col = ind_col)
    
    # Convert FBM to standard matrix and normalize by number of active markers
    K <- K[] / length(ind_col)

    # --- 5. Optimize Lambda (Cross-Validation) ---
    if (is.character(lambda) && lambda == "auto") {
        message("Optimizing Lambda via ", k_folds, "-Fold CV...")

        cand_lambdas <- 10^seq(-3, 3, by = 0.5)
        folds <- cut(seq(1, n), breaks = k_folds, labels = FALSE)

        # Handle Parallel Plan
        if (n_cores > 1) {
            plan(multisession, workers = n_cores)
        } else {
            plan(sequential)
        }

        cv_errors <- future_lapply(cand_lambdas, function(lam) {
            err_sum <- 0
            for (k in 1:k_folds) {
                idx_val <- which(folds == k)
                idx_trn <- which(folds != k)

                K_trn <- K[idx_trn, idx_trn]
                K_val_trn <- K[idx_val, idx_trn]
                y_trn <- y_centered[idx_trn]
                y_val <- y_centered[idx_val]

                # Ridge Solution: alpha = (K + lambda*I)^-1 y
                # Note: We add lambda to diagonal. 
                # Since K is normalized by P, lambda acts as noise/genetic variance ratio.
                I_trn <- diag(length(idx_trn))
                alpha <- solve(K_trn + lam * I_trn, y_trn)

                y_pred <- K_val_trn %*% alpha
                err_sum <- err_sum + sum((y_val - y_pred)^2)
            }
            return(err_sum / n)
        }, future.seed = TRUE)

        # Reset plan
        plan(sequential)

        best_idx <- which.min(unlist(cv_errors))
        lambda <- cand_lambdas[best_idx]
        message("Optimal Lambda selected: ", lambda)
    }

    # --- 6. Final Model Fit ---
    message("Fitting final model...")
    I <- diag(n)
    alpha <- solve(K + lambda * I, y_centered)

    # --- 7. Back-solve for SNP Effects ---
    # u = Z' alpha
    # Since Z was scaled, we must reverse the scaling algebra:
    # u_hat = (raw_prod - center * sum(alpha)) / scale
    
    message("Back-calculating marker effects...")

    # Retrieve scaling parameters used for Z
    scaling <- bigstatsr::big_scale()(geno_mat, ind.col = ind_col)
    centers <- scaling$center
    scales <- scaling$scale

    # raw_prod = X' %*% alpha
    raw_prod <- bigstatsr::big_cprodVec(geno_mat, alpha, ind.col = ind_col)

    sum_alpha <- sum(alpha)
    u_hat <- (raw_prod - centers * sum_alpha) / scales

    # Map effects back to full marker vector (fill monomorphic with 0)
    full_effects <- numeric(p)
    full_effects[ind_col] <- u_hat

    results <- list(
        marker_effects = full_effects,
        model_info = list(
            lambda = lambda,
            intercept = y_mean,   # The mean we subtracted earlier
            n_markers = length(ind_col)
        )
    )
    
    message("Estimation complete.")
    return(results)
}
