# Script: 04_local_gebv.R
# Description: Calculates Local Genomic Estimated Breeding Values (GEBVs) for defined blocks.
# Dependencies: bigstatsr, future.apply

library(bigstatsr)
library(future.apply)

#' Calculate Local GEBVs
#'
#' @param geno_mat Genotype Matrix or FBM.
#' @param marker_effects Vector of marker effects (length p).
#' @param blocks_df Data frame of blocks (BlockID, Start, End).
#' @param n_cores Number of cores.
#' @return A list containing:
#'   - matrix: Local GEBV matrix (Individuals x Blocks).
#'   - variances: Vector of variances per block.
calculate_local_gebv_func <- function(geno_mat, marker_effects, blocks_df, n_cores = 1) {
    n_blocks <- nrow(blocks_df)
    message("Calculating Local GEBVs for ", n_blocks, " blocks...")

    # Pre-calculate scaling stats globally to avoid re-computing inside loop
    # This is much faster
    message("Pre-calculating global scaling statistics...")
    scaling <- bigstatsr::big_scale()(geno_mat)
    centers_all <- scaling$center
    scales_all <- scaling$scale

    if (n_cores > 1) {
        future::plan(future::multisession, workers = n_cores)
    }

    results <- future.apply::future_lapply(1:n_blocks, function(i) {
        start <- blocks_df$Start[i]
        end <- blocks_df$End[i]
        indices <- start:end

        # Extract relevant params
        u_block <- marker_effects[indices]
        mus <- centers_all[indices]
        sigmas <- scales_all[indices]

        # Calculate weights for linear combination
        # GEBV = Z * u
        # Z = (X - mu) / sigma
        # GEBV = (X * u / sigma) - sum(mu * u / sigma)

        weights <- u_block / sigmas
        offset <- sum(mus * weights)

        # Compute X * weights efficiently
        # big_prodVec computes X %*% weights
        gebv_vec <- bigstatsr::big_prodVec(geno_mat, weights, ind.col = indices)
        gebv_vec <- gebv_vec - offset

        return(list(gebv = gebv_vec, var = var(gebv_vec)))
    }, future.seed = TRUE)

    # Aggregate
    local_gebvs <- do.call(cbind, lapply(results, function(x) x$gebv))
    colnames(local_gebvs) <- paste0("Block_", 1:n_blocks)

    variances <- unlist(lapply(results, function(x) x$var))

    return(list(
        matrix = local_gebvs,
        variances = variances
    ))
}
