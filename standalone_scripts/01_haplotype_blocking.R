# Script: 01_haplotype_blocking.R
# Description: Defines haplotype blocks based on Linkage Disequilibrium (LD) decay.
# Dependencies: bigstatsr, data.table

library(data.table)
library(bigstatsr)

#' Define Haplotype Blocks using LD Scanning
#'
#' @param geno_mat A genotype matrix or Filebacked Big Matrix (FBM).
#' @param map_data Optional data.frame containing map information (must have 'chr' column).
#' @param r2_threshold LD r^2 threshold to define block boundaries (default: 0.5).
#' @param window_size Maximum search window size in markers (default: 2000).
#' @param tolerance Number of consecutive SNPs below r2_threshold allowed before breaking LD block (default: 2).
#' @param verbose Logical, print progress messages.
#'
#' @return A data.frame with columns BlockID, Start, End.
define_blocks_ld <- function(geno_mat, map_data = NULL, r2_threshold = 0.5, window_size = 2000, tolerance = 2, verbose = TRUE) {
    if (verbose) message("Starting Haplotype Block Definition (Method: LD)...")

    n_markers <- ncol(geno_mat)
    
    # 1. Handle Chromosome Boundaries
    chr_breaks <- integer(0)
    if (!is.null(map_data) && "chr" %in% names(map_data)) {
        chrs <- map_data$chr
        # Find indices where chr changes (Start of new chromosome)
        if (length(unique(chrs)) > 1) {
            chr_breaks <- which(diff(as.numeric(as.factor(chrs))) != 0) + 1
        }
    }

    starts <- integer()
    ends <- integer()
    current_idx <- 1

    if (verbose) prog_bar <- txtProgressBar(min = 0, max = n_markers, style = 3)

    # 2. Scanning Loop
    while (current_idx <= n_markers) {
        # Define search horizon
        limit_idx <- min(current_idx + window_size, n_markers)

        # Check for chromosome break
        # We find the *next* break that occurs strictly after current_idx but within the limit
        next_break <- chr_breaks[chr_breaks > current_idx & chr_breaks <= limit_idx]
        if (length(next_break) > 0) {
            limit_idx <- min(next_break) - 1
        }

        # Edge case: If limit < current (e.g., immediate break), close block and move on
        if (limit_idx < current_idx) {
            starts <- c(starts, current_idx)
            ends <- c(ends, current_idx)
            current_idx <- current_idx + 1
            if (verbose) setTxtProgressBar(prog_bar, current_idx)
            next
        }

        # Extract window
        window_indices <- current_idx:limit_idx

        # If window is a single marker, it's automatically a block
        if (length(window_indices) < 2) {
            starts <- c(starts, current_idx)
            ends <- c(ends, limit_idx)
            current_idx <- limit_idx + 1
            if (verbose) setTxtProgressBar(prog_bar, current_idx)
            next
        }

        # Extract chunk to memory (Handle FBM vs Matrix)
        if (inherits(geno_mat, "FBM")) {
            mat <- geno_mat[, window_indices, drop = FALSE]
        } else {
            mat <- geno_mat[, window_indices, drop = FALSE]
        }

        # Calculate LD (r^2) against the "Anchor" SNP (first in window)
        r_vals <- cor(mat[, 1], mat, use = "pairwise.complete.obs")
        r2_vals <- as.vector(r_vals^2)
        
        # Handle NAs (e.g. invariant columns) by treating them as 0 correlation
        r2_vals[is.na(r2_vals)] <- 0

        # Scan for block break
        failures <- 0
        block_len <- 0 # 0 means only the anchor itself is kept so far

        # Iterate starting from 2nd marker in window (1st is anchor, r2=1)
        for (j in 2:length(r2_vals)) {
            val <- r2_vals[j]
            
            if (val >= r2_threshold) {
                failures <- 0
                block_len <- j # j is relative index (1-based)
            } else {
                failures <- failures + 1
            }

            if (failures > tolerance) break
        }

        # If block_len is still 0 (immediate failure), the block is just the anchor (length 1)
        if (block_len == 0) block_len <- 1

        final_end <- current_idx + block_len - 1

        # Safety check
        if (final_end < current_idx) final_end <- current_idx

        # Append results
        starts <- c(starts, current_idx)
        ends <- c(ends, final_end)

        # Advance
        current_idx <- final_end + 1
        if (verbose) setTxtProgressBar(prog_bar, current_idx)
    }

    if (verbose) close(prog_bar)

    blocks_df <- data.frame(
        BlockID = seq_along(starts),
        Start = starts,
        End = ends
    )

    if (verbose) message("\nDefined ", nrow(blocks_df), " blocks.")

    return(blocks_df)
}
