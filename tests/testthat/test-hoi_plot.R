test_that("plot_hoi_distribution works", {
    # Initialize object
    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)

    # Mock necessary data
    n_ind <- 20
    n_blocks <- 5

    # Haplo mat (random 1s and 2s)
    haplo$haplo_mat <- matrix(sample(1:2, n_ind * n_blocks, replace = TRUE), nrow = n_ind, ncol = n_blocks)

    # Local GEBVs (Effect distributions)
    local_gebv <- matrix(rnorm(n_ind * n_blocks), nrow = n_ind, ncol = n_blocks)
    haplo$local_gebv <- list(matrix = local_gebv)

    # Test basic run
    expect_error(res <- haplo$plot_hoi_distribution(block_id = 1), NA)

    # Check returns invisible analysis result
    expect_true(is.list(res))
    expect_true(all(c("peaks", "nadir", "p_value", "hoi_haplotypes", "stats") %in% names(res)))

    # Test with checks (indices)
    expect_error(haplo$plot_hoi_distribution(block_id = 1, highlight_ind = c(1, 5)), NA)

    # Test with checks (names with rownames present)
    rownames(local_gebv) <- paste0("Ind", 1:n_ind)
    haplo$local_gebv <- list(matrix = local_gebv)
    expect_error(haplo$plot_hoi_distribution(block_id = 1, highlight_ind = c("Ind1", "Ind5")), NA)

    # Test with checks (names with NO rownames - warning expected)
    haplo$local_gebv <- list(matrix = matrix(rnorm(n_ind * n_blocks), nrow = n_ind)) # Reset without names
    expect_warning(haplo$plot_hoi_distribution(block_id = 1, highlight_ind = c("Ind1")), "Could not match character names")

    # Test validation
    expect_error(haplo$plot_hoi_distribution(block_id = 999), "Invalid block_id")

    haplo$local_gebv <- NULL
    expect_error(haplo$plot_hoi_distribution(block_id = 1), "Local GEBVs not calculated")
})
