test_that("Statistical Engine works (PVE, FA, Model Corr)", {
    # 1. Setup Data
    n_ind <- 50
    n_markers <- 20

    # Random Genotypes
    geno <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)
    rownames(geno) <- paste0("S", 1:n_ind)

    # Phenotype (some correlation with marker 1 and 2)
    pheno <- geno[, 1] * 2 + geno[, 2] + rnorm(n_ind)

    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(geno)
    haplo$load_pheno(pheno)

    # Map
    map <- data.frame(chr = 1, id = paste0("M", 1:n_markers), pos = 1:n_markers, ref = "A", alt = "B")
    haplo$load_map(map)

    # Blocks & Effects
    haplo$define_haploblocks(method = "fixed", window_size = 2) # 10 blocks
    # Mock effects to skip estimation
    haplo$marker_effects <- runif(n_markers)

    # Local GEBV
    haplo$calculate_local_gebv()

    # Signif (mock P_Values)
    haplo$test_significance()

    # 2. Test calculate_pve
    expect_error(haplo$calculate_pve(), NA)
    expect_true(!is.null(haplo$significance$PVE))
    expect_true(all(haplo$significance$PVE >= 0 & haplo$significance$PVE <= 1))

    # 3. Test analyze_block_structure
    expect_error(haplo$analyze_block_structure(top_n = 5, factors = 2), NA)

    # Access private field via hack or add accessor?
    # R6 private fields are locked.
    # We should probably expose `get_fa_results` or similar if we want users to see it for plotting.
    # But the plotting functions are internal/methods of the class, so they can see it.

    # For testing, we can check if print method mentions it or add a temporary accessor.
    # Or rely on `get_model_correlation` which uses it.

    # 4. Test get_model_correlation
    G_smooth <- haplo$get_model_correlation()
    expect_true(is.matrix(G_smooth))
    expect_equal(dim(G_smooth), c(5, 5)) # top_n=5

    # 5. Test analyze_structural_drivers
    res <- haplo$analyze_structural_drivers()
    expect_s3_class(res, "summary.lm")
})
