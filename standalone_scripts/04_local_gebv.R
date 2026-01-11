# Script: 04_local_gebv.R
# Description: Calculates Local Genomic Estimated Breeding Values (GEBVs).
#              Formula: GEBV = Z * u = (X - center)/scale * u
# Dependencies: bigstatsr, future, future.apply

library(bigstatsr)
library(future)
library(future.apply)

#' Calculate Local GEBVs
#'
#' @param geno_mat Genotype Matrix or FBM.
#' @param marker_effects Vector of marker effects (length p).
#' @param blocks_df Data frame of blocks (BlockID, Start, End).
#' @param n_cores Number of cores.
#' @return A list containing:
#'    - matrix: Local GEBV matrix (Individuals x Blocks).
#'    - variances: Vector of variances per block.
calculate_local_gebv_func <- function(geno_mat, marker_effects, blocks_df, n_cores = 1) {
    
    # --- 1. Robust Input Handling ---
    if (!inherits(geno_mat, "FBM")) {
        message("Converting genotype matrix to FBM format...")
        geno_mat <- bigstatsr::as_FBM(geno_mat)
    }

    n_blocks <- nrow(blocks_df)
    n_markers <- ncol(geno_mat)

    # Sanity Check
    if (length(marker_effects) != n_markers) {
        stop("Error: Length of marker_effects (", length(marker_effects), 
             ") does not match number of markers (", n_markers, ").")
    }

    message("Calculating Local GEBVs for ", n_blocks, " blocks (Cores: ", n_cores, ")...")

    # --- 2. Pre-calculate Scaling Statistics ---
    # We scan the whole matrix once to get means (centers) and SDs (scales)
    message("Pre-calculating global scaling statistics...")
    stats <- bigstatsr::big_scale()(geno_mat)
    centers_all <- stats$center
    scales_all <- stats$scale

    # Setup Parallel Plan
    if (n_cores > 1) {
        plan(multisession, workers = n_cores)
    } else {
        plan(sequential)
    }

    # --- 3. Parallel Block Calculation ---
    results <- future_lapply(1:n_blocks, function(i) {
        start <- blocks_df$Start[i]
        end <- blocks_df$End[i]
        indices <- start:end

        # Extract parameters for this block
        u_block <- marker_effects[indices]
        mus <- centers_all[indices]
        sigmas <- scales_all[indices]

        # --- SAFETY CHECK: Monomorphic Markers ---
        # If sigma is ~0, division causes Inf/NaN.
        # Markers with 0 variance typically have 0 effect (u=0), but 0/0 is NaN.
        # We replace invalid sigmas with 1 (or any non-zero) and ensure u is 0.
        bad_idx <- which(sigmas < 1e-8)
        if (length(bad_idx) > 0) {
            sigmas[bad_idx] <- 1.0 
            u_block[bad_idx] <- 0.0 # Force effect to 0 just in case
        }

        # Optimized GEBV Calculation:
        # GEBV = Z u = ((X - mu) / sigma) * u
        #      = (X * u/sigma) - sum(mu * u/sigma)
        
        weights <- u_block / sigmas
        offset <- sum(mus * weights)

        # big_prodVec computes X %*% weights very efficiently
        gebv_vec <- bigstatsr::big_prodVec(geno_mat, weights, ind.col = indices)
        gebv_vec <- gebv_vec - offset

        return(list(gebv = gebv_vec, var = var(gebv_vec)))
    }, future.seed = TRUE)

    # Reset Plan
    plan(sequential)

    # --- 4. Aggregation ---
    # Extract columns
    local_gebvs <- do.call(cbind, lapply(results, function(x) x$gebv))
    colnames(local_gebvs) <- paste0("B", blocks_df$BlockID)

    variances <- unlist(lapply(results, function(x) x$var))

    message("Local GEBV calculation complete.")

    return(list(
        matrix = local_gebvs,
        variances = variances
    ))
}
