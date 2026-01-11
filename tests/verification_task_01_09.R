library(testthat)
library(data.table)
library(bigstatsr)

# Source scripts
source("standalone_scripts/01_haplotype_blocking.R")
source("standalone_scripts/09_vis_association.R")

context("Verification Task: Blocking & Vis")

test_that("define_blocks_ld respects min_physical_dist", {
    # Create Mock Data
    # 10 Markers
    n_mk <- 10
    n_ind <- 20

    # Genotype:
    # Block 1 starts at 1. We want LD to drop at 3.
    # Marker 1 & 2 perfectly correlated.
    # Marker 3 & 4 perfectly correlated but UNCORRELATED to 1 & 2.
    # So LD break should be between 2 and 3.

    X <- matrix(0, nrow = n_ind, ncol = n_mk)
    # Group A (1, 2)
    X[, 1] <- sample(0:2, n_ind, replace = TRUE)
    X[, 2] <- X[, 1] # Perfect LD

    # Group B (3, 4) - independent of A
    X[, 3] <- sample(0:2, n_ind, replace = TRUE)
    X[, 4] <- X[, 3]

    # Map Data
    # Case A: 2 and 3 are Far apart (Regular LD break should happen)
    map_far <- data.frame(chr = rep(1, n_mk), pos = seq(100, by = 10000, length.out = n_mk))
    # Dist(2,3) = 20100 - 10100 = 10000.

    # Case B: 2 and 3 are Close (Should FORCE merge)
    map_close <- data.frame(chr = rep(1, n_mk), pos = seq(100, by = 10000, length.out = n_mk))
    map_close$pos[3] <- map_close$pos[2] + 10 # Marker 3 is only 10bp after Marker 2

    # Test 1: Standard LD blocking (High threshold -> should break at 3)
    cat("\nRunning Standard LD Block...\n")
    blocks_std <- define_blocks_ld(X, map_far, r2_threshold = 0.9, min_physical_dist = 0, verbose = FALSE)

    # Expect block 1 to end at 2 (since 2 and 3 are uncorrelated)
    expect_equal(blocks_std$End[1], 2)
    expect_equal(blocks_std$Start[2], 3)

    # Test 2: Physical Constraint (Minimum Block Size)
    # We want to verify that blocks are forced to continue if length < min_physical_dist.

    # Setup Map C:
    # 1: 100
    # 2: 105 (R2 Low w.r.t 1, but Dist=5 < 50) -> Should Keep
    # 3: 110 (R2 Low w.r.t 1, but Dist=10 < 50) -> Should Keep
    # 4: 200 (R2 Low w.r.t 1, Dist=100 > 50) -> Should Break

    X_min <- matrix(0, nrow = n_ind, ncol = 4)
    # All independent (R2 low)
    X_min[, 1] <- sample(0:2, n_ind, replace = TRUE)
    X_min[, 2] <- sample(0:2, n_ind, replace = TRUE)
    X_min[, 3] <- sample(0:2, n_ind, replace = TRUE)
    X_min[, 4] <- sample(0:2, n_ind, replace = TRUE)

    # Ensure variance to avoid cor warnings
    for (k in 1:4) {
        if (sd(X_min[, k]) == 0) X_min[1, k] <- (X_min[1, k] + 1) %% 3
        if (sd(X_min[, k]) == 0) X_min[2, k] <- (X_min[2, k] + 1) %% 3
    }

    map_min <- data.frame(chr = rep(1, 4), pos = c(100, 105, 110, 200))

    cat("\nRunning Min Dist Block...\n")
    blocks_min <- define_blocks_ld(X_min, map_min, r2_threshold = 0.9, min_physical_dist = 50, verbose = FALSE)

    # Block 1 should start at 1.
    # Marker 2: Dist 5. Keep.
    # Marker 3: Dist 10. Keep.
    # Marker 4: Dist 100. Break.
    # So Block 1 should be 1-3. End = 3.

    expect_equal(blocks_min$End[1], 3)
})

test_that("plot_haplo_freq_trend runs correctly", {
    # Setup Mock Data
    sup <- data.frame(BlockID = 1, BestHaploID = 1, EffectSize = 0.5)
    haplo <- matrix(c(1, 1, 1, 0, 0, 0), ncol = 1) # 6 indivs
    rownames(haplo) <- paste0("Ind", 1:6)
    meta <- data.frame(ID = paste0("Ind", 1:6), Year = rep(2020:2022, each = 2))

    # Should run without error
    pdf(file = NULL)
    expect_silent(plot_haplo_freq_trend(sup, haplo, meta, top_k = 1))
    dev.off()
})
