test_that("scale_haplo_effects centers columns to mean 0", {
    # Initialize object
    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)

    # Create a dummy matrix
    set.seed(123)
    local_gebv <- matrix(rnorm(100, mean = 5, sd = 1), ncol = 10)

    # Test with explicit argument
    scaled_gebv <- haplo$scale_haplo_effects(local_gebv)

    # Check dimensions
    expect_equal(dim(scaled_gebv), dim(local_gebv))

    # Check column means (should be close to 0)
    col_means <- colMeans(scaled_gebv)
    expect_true(all(abs(col_means) < 1e-10))

    # Check if relative differences are preserved
    expect_equal(sd(scaled_gebv[, 1]), sd(local_gebv[, 1]))

    # Test with internal data
    haplo$local_gebv <- list(matrix = local_gebv) # Mocking internal state
    scaled_internal <- haplo$scale_haplo_effects() # No args
    expect_equal(scaled_internal, scaled_gebv)
})

test_that("analyze_hoi detects peaks and superior haplotypes", {
    # Initialize object
    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)

    n_ind <- 100
    n_blocks <- 1

    # Create haplo_geno: 1s and 2s
    haplo_geno <- matrix(sample(1:2, n_ind, replace = TRUE), ncol = 1)
    haplo$haplo_mat <- haplo_geno # Mocking internal state

    # Create local_gebv based on haplotypes
    # H1 -> effect around -2, H2 -> effect around 2
    effects <- numeric(n_ind)
    effects[haplo_geno == 1] <- rnorm(sum(haplo_geno == 1), mean = -2, sd = 0.5)
    effects[haplo_geno == 2] <- rnorm(sum(haplo_geno == 2), mean = 2, sd = 0.5)

    local_gebv <- matrix(effects, ncol = 1)
    haplo$local_gebv <- list(matrix = local_gebv) # Mocking internal state

    # Run analysis with explicit args
    res <- haplo$analyze_hoi(block_id = 1, local_gebv = local_gebv)

    # Checks
    expect_type(res, "list")
    expect_named(res, c("peaks", "nadir", "p_value", "hoi_haplotypes", "stats"))

    # We expect 2 peaks
    expect_length(res$peaks, 2)
    expect_true(res$peaks[1] < 0) # Low peak
    expect_true(res$peaks[2] > 0) # High peak

    # p-value should be significant
    expect_true(res$p_value < 0.05)

    # HOI haplotypes should include 2 (the high effect one)
    expect_true(2 %in% res$hoi_haplotypes)
    expect_false(1 %in% res$hoi_haplotypes)

    # Run analysis with internal data
    res_internal <- haplo$analyze_hoi(block_id = 1)
    expect_equal(res_internal, res)

    # Warning if too few individuals/peaks (e.g. unimodal)
    local_gebv_uni <- matrix(rnorm(100), ncol = 1)
    expect_warning(haplo$analyze_hoi(block_id = 1, local_gebv = local_gebv_uni), "Less than 2 peaks detected")
})
