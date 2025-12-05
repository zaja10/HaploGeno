test_that("scale_haplo_effects centers columns to mean 0", {
    # Create a dummy matrix
    set.seed(123)
    local_gebv <- matrix(rnorm(100, mean = 5, sd = 1), ncol = 10)

    scaled_gebv <- scale_haplo_effects(local_gebv)

    # Check dimensions
    expect_equal(dim(scaled_gebv), dim(local_gebv))

    # Check column means (should be close to 0)
    col_means <- colMeans(scaled_gebv)
    expect_true(all(abs(col_means) < 1e-10))

    # Check if relative differences are preserved (variance shouldn't change, just shifted)
    expect_equal(sd(scaled_gebv[, 1]), sd(local_gebv[, 1]))
})

test_that("analyze_hoi detects peaks and superior haplotypes", {
    # Mock HaploObject
    # We need to simulate haplo_geno matrix
    # Let's say we have 100 individuals, and 1 block
    # We'll simulate 2 haplotypes: H1 (low effect) and H2 (high effect)

    n_ind <- 100
    n_blocks <- 1

    # Create haplo_geno: 1s and 2s
    haplo_geno <- matrix(sample(1:2, n_ind, replace = TRUE), ncol = 1)

    haplo_obj <- list(haplo_geno = haplo_geno)

    # Create local_gebv based on haplotypes
    # H1 -> effect around -2, H2 -> effect around 2
    effects <- numeric(n_ind)
    effects[haplo_geno == 1] <- rnorm(sum(haplo_geno == 1), mean = -2, sd = 0.5)
    effects[haplo_geno == 2] <- rnorm(sum(haplo_geno == 2), mean = 2, sd = 0.5)

    local_gebv <- matrix(effects, ncol = 1)

    # Run analysis
    res <- analyze_hoi(haplo_obj, local_gebv, block_id = 1)

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

    # Check if it handles unimodal case gracefully (or just returns warning/simple stats)
    # This depends on implementation details, but let's test a simple uniform case
    # where it might not find 2 clear peaks or warnings.
    # Our implementation warns if < 2 peaks.

    local_gebv_uni <- matrix(rnorm(100), ncol = 1)
    expect_warning(res_uni <- analyze_hoi(haplo_obj, local_gebv_uni, block_id = 1), "Less than 2 peaks detected")
    expect_null(res_uni$hoi_haplotypes)
})
