# Script: 05_significance_testing.R
# Description: Assesses significance of local genomic regions using Empirical Permutation.
#              Generates a null distribution of local variances by permuting phenotypes.
# Dependencies: bigstatsr, stats

library(bigstatsr)
library(stats)

#' Test Significance of Local Blocks
#'
#' @param geno_mat Genotype Matrix or FBM.
#' @param pheno_vec Phenotype vector.
#' @param local_gebv_list Result from calculate_local_gebv_func (must have 'variances').
#' @param blocks_df Data frame of blocks.
#' @param n_perms Number of permutations (default: 10).
#' @param n_null_blocks Number of random blocks to sample per permutation (default: 50).
#' @param lambda Optional lambda used in model training (default: 1.0).
#' @return Data frame with BlockID, Variance, P_Value.
test_significance_func <- function(geno_mat, pheno_vec, local_gebv_list, blocks_df,
                                   n_perms = 10, n_null_blocks = 50, lambda = 1.0) {
    
    # --- 1. Robust Input Handling ---
    if (!inherits(geno_mat, "FBM")) {
        message("Converting genotype matrix to FBM format...")
        geno_mat <- bigstatsr::as_FBM(geno_mat)
    }

    n <- nrow(geno_mat)
    
    # Center phenotypes to ensure permutations don't struggle with intercepts
    y_centered <- as.vector(pheno_vec) - mean(pheno_vec, na.rm = TRUE)

    message("Running Significance Test (Empirical Permutation)...")

    # --- 2. Setup Null Distribution ---
    # Identify polymorphic markers for consistent K calculation
    stats <- bigstatsr::big_colstats(geno_mat)
    ind_col <- which(stats$var > 1e-8)
    n_active <- length(ind_col)

    # Pre-calc Scaling (Centers/Scales)
    message("Pre-calculating global scaling statistics...")
    scaling <- bigstatsr::big_scale()(geno_mat)
    centers_all <- scaling$center
    scales_all <- scaling$scale

    # Pre-calc Kernel K for fast solving inside permutation loop
    message("Computing and Inverting Kernel (N=", n, ") for fast permutation...")
    
    # Calculate K using only polymorphic markers
    K <- bigstatsr::big_tcrossprodSelf(geno_mat, fun.scaling = bigstatsr::big_scale(), ind.col = ind_col)
    K <- K[] / n_active

    # Inverse of (K + lambda I)
    # We invert this ONCE. For each permutation, we just multiply Ki %*% y_perm
    Ki <- solve(K + lambda * diag(n))

    null_vars <- numeric(n_perms * n_null_blocks)
    counter <- 1

    # --- 3. Permutation Loop ---
    # Note: This loop is typically fast enough to run sequentially in R 
    # because the heavy lifting (Ki) is already done.
    
    for (i in seq_len(n_perms)) {
        if (i %% 5 == 0) message(sprintf("  Permutation %d/%d...", i, n_perms))
        
        # Shuffle Phenotypes (break genetic link)
        y_perm <- sample(y_centered)

        # Fast Alpha Calculation: alpha = Ki * y_perm
        alpha_perm <- Ki %*% y_perm
        sum_alpha <- sum(alpha_perm)

        # Sample random blocks to estimate null variance
        # We don't need to calculate every block for every permutation (too slow)
        rand_blk_ids <- sample(seq_len(nrow(blocks_df)), min(n_null_blocks, nrow(blocks_df)))

        for (b_id in rand_blk_ids) {
            start <- blocks_df$Start[b_id]
            end <- blocks_df$End[b_id]
            indices <- start:end

            # Safety: Handle monomorphic markers in specific blocks
            # If scale is ~0, avoid division by zero
            blk_scales <- scales_all[indices]
            bad_idx <- which(blk_scales < 1e-8)
            if (length(bad_idx) > 0) blk_scales[bad_idx] <- 1.0

            # Estimate null marker effects: u = Z' alpha
            # u = (X' alpha - sum(alpha)*mu) / sigma
            
            # X' alpha for specific block columns
            raw_prod <- bigstatsr::big_cprodVec(geno_mat, alpha_perm, ind.col = indices)
            u_hat <- (raw_prod - centers_all[indices] * sum_alpha) / blk_scales
            
            # Force 0 effect for monomorphic markers
            if (length(bad_idx) > 0) u_hat[bad_idx] <- 0

            # Calculate GEBV for this null effect
            # GEBV = Z * u = (X - mu)/sigma * u
            weights <- u_hat / blk_scales
            offset <- sum(centers_all[indices] * weights)
            
            gebv_null <- bigstatsr::big_prodVec(geno_mat, weights, ind.col = indices) - offset

            null_vars[counter] <- var(gebv_null)
            counter <- counter + 1
        }
    }

    # --- 4. Calculate P-values ---
    obs_vars <- local_gebv_list$variances

    # ECDF of null distribution
    # p = 1 - ECDF(obs)
    # This represents the probability of observing a variance >= obs under the null
    null_ecdf <- ecdf(null_vars)
    p_values <- 1 - null_ecdf(obs_vars)

    # Avoid exactly 0 p-values (for -log10 plots)
    # Set 0s to 1 / (total_null_samples + 1)
    min_p <- 1 / (length(null_vars) + 1)
    p_values[p_values == 0] <- min_p

    results <- data.frame(
        BlockID = blocks_df$BlockID,
        Variance = obs_vars,
        P_Value = p_values
    )

    message("Significance testing complete.")
    return(results)
}
