test_that("Selection and Profiling Tools run without error", {
    # Setup Mock Data
    n_ind <- 20
    n_markers <- 100
    geno <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)
    rownames(geno) <- paste0("S", 1:n_ind)

    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(geno)

    # Map
    map <- data.frame(chr = 1, id = paste0("M", 1:n_markers), pos = 1:n_markers, ref = "A", alt = "B")
    haplo$load_map(map)

    # Pheno
    # Create a dummy phenotype correlated with first few markers to ensure high variance
    eff <- rep(0, n_markers)
    eff[1:10] <- 1
    pheno <- (geno %*% eff) + rnorm(n_ind)
    haplo$load_pheno(as.vector(pheno))

    # Blocks & Effects
    haplo$define_haploblocks(method="fixed", window_size = 5)
    haplo$marker_effects <- runif(n_markers)

    # Compute Required Stats
    haplo$encode_haplotypes()
    haplo$calculate_local_gebv()
    haplo$test_significance() # Required for selection

    # 1. Test get_haplo_extremes
    extremes <- haplo$get_haplo_extremes(top_n = 5, threshold = 0.5)
    expect_true(is.list(extremes))
    expect_equal(length(extremes), 5)

    # Check structure of one result
    res1 <- extremes[[1]]
    expect_true(!is.null(res1$Positive))
    expect_true(!is.null(res1$Negative))
    expect_true(is.numeric(res1$Mean))

    # 2. Test plot_haplo_profile
    pdf(NULL)

    # Standard run
    expect_error(haplo$plot_haplo_profile(top_n_blocks = 5, sort_by = "pheno"), NA)

    # With subset
    subset_ids <- rownames(geno)[1:10]
    expect_error(haplo$plot_haplo_profile(genotypes = subset_ids, top_n_blocks = 5), NA)

    # With numeric subset
    expect_error(haplo$plot_haplo_profile(genotypes = 1:5, top_n_blocks = 5), NA)

    dev.off()
})
