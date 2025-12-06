test_that("Active markers bug in calculate_local_gebv is fixed", {
    # 1. Setup small genotype matrix (10 samples, 20 markers)
    n_ind <- 10
    n_markers <- 20
    geno <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)
    rownames(geno) <- paste0("Sample_", 1:n_ind)

    # Introduce some monomorphic markers (e.g. cols 5, 10, 15)
    geno[, 5] <- 0
    geno[, 10] <- 2

    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)

    # 2. Import genotypes
    haplo$import_genotypes(geno)

    # Check if sample_ids are captured (New Requirement)
    # This part of the test might fail initially or verify the need for the feature
    # expect_equal(haplo$sample_ids, rownames(geno))

    # 3. Simulate filter_monomorphic (manual setting of active_markers)
    # Suppose we keep all except 5 and 10.
    # Indices: 1, 2, 3, 4, 6, 7, 8, 9, 11...
    keep_cols <- setdiff(1:n_markers, c(5, 10))
    haplo$active_markers <- keep_cols

    # 4. Load Mock Map (needed for blocks)
    map <- data.frame(
        chr = rep(1, n_markers),
        id = paste0("M", 1:n_markers),
        pos = 1:n_markers,
        ref = "A",
        alt = "T"
    )
    haplo$load_map(map)

    # 5. Define blocks (Fixed window of 5 markers)
    # But 5 markers in REMOTE/RELATIVE index space (of 18 active markers)
    # Blocks should be based on length(active_markers) = 18.
    haplo$define_blocks_fixed(window_size = 5)

    # 6. Mock marker effects (size of active_markers)
    haplo$marker_effects <- rnorm(length(keep_cols))

    # 7. Calculate local GEBV
    # This should NOT fail and should produce non-zero variance results if data is random enough
    # If the bug exists, it might assume block indices (1:5) match columns 1:5, but if we dropped col 1, it would mismatch.
    # To make it tricky, let's drop the FIRST column.

    # RETRY with specific scenario: Disjoint indices.
    # Map cols 11-20 to active markers.
    keep_cols_3 <- 11:n_markers
    haplo$active_markers <- keep_cols_3
    haplo$marker_effects <- rep(1, length(keep_cols_3)) # All effects = 1
    haplo$define_blocks_fixed(window_size = 5)

    # Block 1: Rel 1-5. Active: 11, 12, 13, 14, 15.
    # Current Buggy Code: Check 1:5 %in% 11:20 -> None. Returns 0.
    # Correct Code: Map 1:5 -> 11:15. Returns sum of effects (approx).

    haplo$calculate_local_gebv()

    # Check Block 1 result
    # access matrix directly. Block 1 is column 1.
    res_block_1 <- haplo$local_gebv$matrix[, 1]

    # If bug exists, this will be all zeros
    # We expect non-zeros because effects are 1 and genotypes are 0/1/2.
    # Unless all genotypes for these markers are 0 (unlikely with random sampling)
    expect_true(sum(abs(res_block_1)) > 1e-5)

    # Check if checking for valid_idx_in_block worked correctly.
    # If bug existed: indices = 1:5. active_markers = 2:20.
    # valid_idx_in_block = which(1:5 %in% 2:20) -> which(c(F, T, T, T, T)) -> 2,3,4,5.
    # abs_indices = indices[2,3,4,5] = 2,3,4,5.
    # real_indices (for big_prod) = abs_indices = 2,3,4,5.
    # BUT we wanted to access the FIRST 5 active markers, which are columns 2,3,4,5,6.
    # The code effectively skipped column 6 (the 5th active marker) and maybe tried to read column 1 (which is active? No col 1 is dropped).
    # Wait, if `indices` are relative (1:5), they are pointers to `active_markers`.
    # `active_markers[1:5]` -> 2,3,4,5,6.
    # The BUGGY code does: `indices %in% active_markers`. 1:5 %in% 2:20.
    # 1 is NOT in 2:20. So it drops index 1.
    # So it computes GEBV using only markers 2,3,4,5.
    # It SHOULD use 2,3,4,5,6.
    # So the result would be wrong (missing one marker effect).

    # We can't easily assert "wrongness" on random data without a ground truth, but we can check matching dimensions
    # and that we processed the expected number of markers.
})

test_that("Sample IDs are persisted", {
    n_ind <- 5
    n_markers <- 10
    geno <- matrix(0, nrow = n_ind, ncol = n_markers)
    rownames(geno) <- LETTERS[1:5]

    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)
    haplo$import_genotypes(geno)

    expect_equal(haplo$sample_ids, LETTERS[1:5])
})

test_that("get_block_markers works", {
    n_ind <- 5
    n_markers <- 20
    geno <- matrix(0, nrow = n_ind, ncol = n_markers)
    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(geno)

    # Load Map
    map <- data.frame(
        chr = rep(1, n_markers),
        id = paste0("M", 1:n_markers),
        pos = 1:n_markers,
        ref = "A",
        alt = "T"
    )
    haplo$load_map(map)

    # Blocks (fixed)
    haplo$define_blocks_fixed(window_size = 5)

    # Check Block 1 markers
    # Default (all active)
    markers_b1 <- haplo$get_block_markers(1)
    expect_equal(markers_b1, paste0("M", 1:5))

    # With active markers subset
    haplo$active_markers <- 11:20
    haplo$define_blocks_fixed(window_size = 5)
    # Block 1 in relative space (1:5) -> Active 11:15 -> Markers M11:M15
    markers_b1_active <- haplo$get_block_markers(1)
    expect_equal(markers_b1_active, paste0("M", 11:15))
})
