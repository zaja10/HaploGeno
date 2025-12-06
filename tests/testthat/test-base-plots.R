test_that("Base R Visualizations run without error", {
    # Setup Mock Data
    n_ind <- 20
    n_markers <- 50
    geno <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)
    rownames(geno) <- paste0("S", 1:n_ind)

    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(geno)

    # Map
    map <- data.frame(chr = 1, id = paste0("M", 1:n_markers), pos = 1:n_markers, ref = "A", alt = "B")
    haplo$load_map(map)

    # Pheno
    haplo$load_pheno(rnorm(n_ind))

    # Blocks & Effects
    haplo$define_blocks_fixed(window_size = 5)
    haplo$marker_effects <- runif(n_markers)

    # Compute HRM & GEBV
    haplo$encode_haplotypes() # Required for HRM
    haplo$compute_hrm()
    haplo$calculate_local_gebv()

    # Run Stats Engine
    haplo$test_significance()
    haplo$analyze_block_structure(top_n = 5, factors = 2)

    # Test Biplot
    pdf(NULL) # Prevent plot window
    expect_error(haplo$plot_haplo_biplot(), NA)

    # Test FA Genome
    expect_error(haplo$plot_fa_genome(), NA)

    # Test Communality
    expect_error(haplo$plot_communality(), NA)

    # Test Scree
    expect_error(haplo$plot_scree(), NA)
    dev.off()
})
