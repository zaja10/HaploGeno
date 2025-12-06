test_that("End-to-end integration with minimal CSV input", {
    # 1. Create minimal CSV (Dosage, no headers, no map)
    # 100 individuals, 200 markers
    n_ind <- 50
    n_markers <- 200

    # Random genotypes 0,1,2
    geno_mat <- matrix(sample(0:2, n_ind * n_markers, replace = TRUE), nrow = n_ind, ncol = n_markers)

    # Phenotypes (Outcome = sum of first 5 markers effects + noise)
    effects <- rep(0, n_markers)
    effects[1:5] <- c(2, -2, 1, -1, 0.5)
    y <- geno_mat %*% effects + rnorm(n_ind)

    # Write to CSV without headers (simulate raw matrix)
    tf_geno <- tempfile(fileext = ".csv")
    write.table(geno_mat, tf_geno, sep = ",", row.names = FALSE, col.names = FALSE)

    # Write pheno
    tf_pheno <- tempfile(fileext = ".txt")
    write.table(y, tf_pheno, row.names = FALSE, col.names = FALSE)

    # 2. Pipeline
    ho <- HaploObject$new(tempfile())

    # Import Genotypes
    # Expect header detection to likely fail or use V1, V2... which is fine.
    # Our robust importer checks if V1 is character. Here it's numeric. So it should read as matrix.
    # It will assume V1, V2 are column names if header=FALSE not specified in fread?
    # fread auto-detects.
    ho$import_genotypes(tf_geno)

    expect_equal(nrow(ho$geno), n_ind)
    expect_equal(ncol(ho$geno), n_markers)

    # Map should be auto-generated
    expect_true(!is.null(ho$map))
    expect_equal(nrow(ho$map), n_markers)

    # Load Pheno
    ho$load_pheno(tf_pheno)
    expect_equal(length(ho$pheno), n_ind)

    # Define Blocks (Fixed)
    ho$define_haploblocks(method = "fixed", window_size = 10)
    expect_true(!is.null(ho$blocks))

    # Estimate Effects
    ho$estimate_marker_effects()
    expect_true(!is.null(ho$marker_effects))

    # Calculate Local GEBV
    ho$calculate_local_gebv()
    expect_true(!is.null(ho$local_gebv))

    # Significance
    ho$test_significance()
    expect_true(!is.null(ho$significance))

    # Plotting (Should not error)
    expect_error(ho$plot_manhattan(), NA)

    # 3. Phased/Character Integration
    # Test with "0/1" format to ensure it flows through pipeline
    tf_phased <- tempfile(fileext = ".csv")
    phased_mat <- matrix("0/0", 10, 10)
    phased_mat[1:5, 1] <- "0/1"
    write.csv(as.data.frame(phased_mat), tf_phased, row.names = FALSE)

    ho2 <- HaploObject$new(tempfile())
    ho2$import_genotypes(tf_phased)
    expect_equal(ho2$geno[1, 1], 1)
    expect_equal(ho2$geno[6, 1], 0)

    # Dummy pheno
    ho2$pheno <- rnorm(10)

    ho2$filter_monomorphic()
    ho2$define_haploblocks(method = "fixed", window_size = 5)
    ho2$estimate_marker_effects()
    ho2$calculate_local_gebv()
    expect_error(ho2$test_significance(), NA)
})
