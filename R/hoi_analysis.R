#' Scale Haplotype Effects
#'
#' Centers the local GEBV matrix columns to have a mean of zero.
#'
#' @param local_gebv A matrix of local GEBVs (N x n_blocks).
#' @return A matrix of scaled local GEBVs.
#' @export
scale_haplo_effects <- function(local_gebv) {
    # Center columns to mean 0
    scaled_gebv <- scale(local_gebv, center = TRUE, scale = FALSE)
    return(scaled_gebv)
}

#' Analyze Haplotypes of Interest (HOI)
#'
#' Identifies superior haplotypes at a specific block by analyzing the distribution
#' of local GEBVs. Assumes a bimodal distribution where one peak represents
#' the HOI.
#'
#' @param haplo_obj A HaploObject with defined blocks and encoded haplotypes.
#' @param local_gebv A matrix of local GEBVs (N x n_blocks).
#' @param block_id The ID of the block to analyze.
#' @return A list containing peak values, p-value, HOI haplotypes, and stats.
#' @importFrom stats density t.test
#' @export
analyze_hoi <- function(haplo_obj, local_gebv, block_id) {
    if (is.null(haplo_obj$haplo_geno)) stop("Haplotypes must be encoded first.")

    # Extract data for the block
    # Check if block_id is valid
    if (!block_id %in% 1:ncol(local_gebv)) stop("Invalid block_id.")

    effects <- local_gebv[, block_id]

    # Scale effects
    scaled_effects <- scale(effects, center = TRUE, scale = FALSE)

    # Estimate density
    d <- stats::density(scaled_effects)

    # Find peaks (local maxima)
    # Simple peak finding: where derivative changes from positive to negative
    # We look for indices where d$y[i] > d$y[i-1] and d$y[i] > d$y[i+1]
    peak_indices <- which(diff(sign(diff(d$y))) == -2) + 1
    peak_x <- d$x[peak_indices]
    peak_y <- d$y[peak_indices]

    if (length(peak_indices) < 2) {
        warning("Less than 2 peaks detected. Returning simple stats.")
        return(list(
            peaks = peak_x,
            p_value = NA,
            hoi_haplotypes = NULL,
            stats = summary(scaled_effects)
        ))
    }

    # Sort peaks by x value to separate "low" and "high" groups
    # We assume the "superior" peak is the one with higher effect value
    sorted_peaks <- sort(peak_x)

    # Identify the two most prominent peaks (highest density)
    # This helps ignore small noisy peaks
    top_peaks_idx <- order(peak_y, decreasing = TRUE)[1:2]
    main_peaks_x <- sort(peak_x[top_peaks_idx])

    low_peak <- main_peaks_x[1]
    high_peak <- main_peaks_x[2]

    # Find nadir between them
    # Range between peaks
    range_indices <- which(d$x > low_peak & d$x < high_peak)
    if (length(range_indices) > 0) {
        nadir_idx <- range_indices[which.min(d$y[range_indices])]
        nadir_x <- d$x[nadir_idx]
    } else {
        nadir_x <- (low_peak + high_peak) / 2 # Fallback
    }

    # Split individuals
    group_high <- which(scaled_effects > nadir_x)
    group_low <- which(scaled_effects <= nadir_x)

    # Significance test
    # t-test between the two groups
    if (length(group_high) > 1 && length(group_low) > 1) {
        test_res <- stats::t.test(scaled_effects[group_high], scaled_effects[group_low], alternative = "greater")
        p_val <- test_res$p.value
    } else {
        p_val <- NA
    }

    # Identify haplotypes associated with the high group
    # We need the haplotype integers for this block
    haplo_col <- haplo_obj$haplo_geno[, block_id]

    # Get unique haplotypes in the high group
    # We aggregate effects by haplotype
    haplo_means <- tapply(scaled_effects, haplo_col, mean)
    hoi_haplos <- as.integer(names(haplo_means[haplo_means > nadir_x]))

    return(list(
        peaks = main_peaks_x,
        nadir = nadir_x,
        p_value = p_val,
        hoi_haplotypes = hoi_haplos,
        stats = list(
            mean_high = mean(scaled_effects[group_high]),
            mean_low = mean(scaled_effects[group_low]),
            n_high = length(group_high),
            n_low = length(group_low)
        )
    ))
}
