# Script: 05_significance_testing.R
# Description: Assesses significance of local genomic regions using Empirical Permutation.
# Dependencies: bigstatsr, stats

library(bigstatsr)

#' Test Significance of Local Blocks
#'
#' @param geno_mat Genotype Matrix/FBM.
#' @param pheno_vec Phenotype vector.
#' @param local_gebv_list Result from calculate_local_gebv_func (has 'variances').
#' @param blocks_df Data frame of blocks.
#' @param n_perms Number of permutations (default: 10).
#' @param n_null_blocks Number of random blocks to sample per permutation (default: 50).
#' @param lambdas Optional lambda used in model training (default: 1.0).
#' @return Data frame with BlockID, Variance, P_Value.
test_significance_func <- function(geno_mat, pheno_vec, local_gebv_list, blocks_df,
                                   n_perms = 10, n_null_blocks = 50, lambda = 1.0) {
    message("Running Significance Test (Empirical Permutation)...")

    # 1. Setup Null Distribution
    null_vars <- numeric(n_perms * n_null_blocks)
    n <- nrow(geno_mat)

    # Pre-calc Scaling
    scaling <- bigstatsr::big_scale()(geno_mat)
    centers_all <- scaling$center
    scales_all <- scaling$scale

    # Pre-calc Kernel K for fast solving inside permutation loop
    message("Computing/Inverting Kernel for fast permutation...")
    K <- bigstatsr::big_tcrossprodSelf(geno_mat, fun.scaling = bigstatsr::big_scale())
    K <- K[] / ncol(geno_mat)

    # Inverse of (K + lambda I)
    # We invert this ONCE so we can quickly calculate alpha for permuted y
    Ki <- solve(K + lambda * diag(n))

    counter <- 1

    # 2. Permutation Loop
    for (i in seq_len(n_perms)) {
        # Shuffle Phenotypes
        y_perm <- sample(pheno_vec)

        # Fast Alpha: alpha = Ki * y_perm
        alpha_perm <- Ki %*% y_perm
        sum_alpha <- sum(alpha_perm)

        # Sample random blocks to estimate null variance
        rand_blk_ids <- sample(seq_len(nrow(blocks_df)), min(n_null_blocks, nrow(blocks_df)))

        for (b_id in rand_blk_ids) {
            start <- blocks_df$Start[b_id]
            end <- blocks_df$End[b_id]
            indices <- start:end

            # Estimate null marker effects for this block: u = Z' alpha
            # u = (X' alpha - sum(alpha)*mu) / sigma (simplified derived math)

            # X' alpha for specific cols
            raw_prod <- bigstatsr::big_cprodVec(geno_mat, alpha_perm, ind.col = indices)
            u_hat <- (raw_prod - centers_all[indices] * sum_alpha) / scales_all[indices]

            # Calculate GEBV for this null effect
            weights <- u_hat / scales_all[indices]
            offset <- sum(centers_all[indices] * weights)
            gebv_null <- bigstatsr::big_prodVec(geno_mat, weights, ind.col = indices) - offset

            null_vars[counter] <- var(gebv_null)
            counter <- counter + 1
        }
    }

    # 3. Calculate P-values
    obs_vars <- local_gebv_list$variances

    # ECDF of null distribution
    null_ecdf <- ecdf(null_vars)
    p_values <- 1 - null_ecdf(obs_vars)

    # Avoid exactly 0 p-values for log plots
    min_p <- 1 / (length(null_vars) + 1)
    p_values[p_values == 0] <- min_p

    results <- data.frame(
        BlockID = blocks_df$BlockID,
        Variance = obs_vars,
        P_Value = p_values
    )

    return(results)
}
