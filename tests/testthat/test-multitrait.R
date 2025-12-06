test_that("Multi-Trait FA works", {
    # 1. Setup Data
    n_ind <- 40
    n_markers <- 20

    # Genotypes
    geno <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)
    rownames(geno) <- paste0("S", 1:n_ind)

    # Phenotypes (Two traits sharing some genetic basis)
    # Factor ~ Marker 1
    f <- geno[, 1]
    y1 <- f * 2 + rnorm(n_ind)
    y2 <- f * -1.5 + rnorm(n_ind) # correlated negatively

    # Object 1
    h1 <- HaploObject$new(tempfile())
    h1$import_genotypes(geno)
    h1$load_pheno(y1)
    h1$define_haploblocks(method = "fixed", window_size = 2)
    h1$marker_effects <- runif(n_markers) # mock
    h1$calculate_local_gebv()

    # Object 2
    h2 <- HaploObject$new(tempfile())
    h2$import_genotypes(geno)
    h2$load_pheno(y2)
    h2$define_haploblocks(method = "fixed", window_size = 2)
    h2$marker_effects <- runif(n_markers) # mock
    h2$calculate_local_gebv()

    # 2. Run Multi-Trait Analysis
    res <- analyze_multitrait_structure(list(h1, h2), factors = 2, top_n = 10)

    # 3. Assertions
    expect_s3_class(res, "MultiTraitFA")
    expect_equal(dim(res$Scores), c(n_ind, 2))
    expect_true(nrow(res$Loadings) > 0)

    # Check trait indices
    expect_true(all(res$TraitIndices %in% c(1, 2)))

    # Print method
    expect_output(print(res), "Multi-Trait Factor Analysis Result")

    # Clean up
    rm(h1, h2)
    gc()
})
