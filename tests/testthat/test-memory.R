test_that("Temporary backing files are cleaned up on object destruction", {
    # Skip if not capable of file operations (e.g. strict CRAN checks on some systems)
    # skip_on_cran() removed for local verification

    # Create a dummy genotype matrix
    n <- 50
    p <- 100
    geno <- matrix(rnorm(n * p), n, p)
    rownames(geno) <- paste0("Sample", 1:n)

    # 1. Define a temporary path
    temp_bk <- tempfile()

    # 2. Instantiate Object
    obj <- HaploObject$new(backing_file_path = temp_bk)

    # 3. Import Genotypes (Should create backing file and mark as temp)
    # using matrix input
    obj$import_genotypes(geno)

    # Verify backing file exists
    full_bk <- paste0(temp_bk, ".bk")
    full_rds <- paste0(temp_bk, ".rds")

    expect_true(file.exists(full_bk))
    # expect_true(file.exists(full_rds)) # as_FBM does not create RDS by default

    # 4. Destroy Object
    rm(obj)
    gc()
    gc() # Force garbage collection (twice)

    # 5. Verify files are gone
    expect_false(file.exists(full_bk))
    # expect_false(file.exists(full_rds))
})

test_that("User-supplied existing files are NOT deleted", {
    # skip_on_cran() removed for local verification

    # Create a "permanent" backing file manually
    n <- 10
    p <- 10
    mat <- matrix(rnorm(n * p), n, p)
    fbm_path <- tempfile()
    FBM <- bigstatsr::as_FBM(mat, backingfile = fbm_path)
    # Save the FBM descriptor to mimic an existing RDS
    rds_path <- paste0(fbm_path, ".rds")
    bk_path <- paste0(fbm_path, ".bk")
    saveRDS(FBM, rds_path)

    expect_true(file.exists(bk_path))

    # Instantiate HaploObject pointing to this existing file
    # We pass the RDS path to import_genotypes
    obj <- HaploObject$new(backing_file_path = tempfile()) # Path here won't be used if we load existing

    # We need to simulate loading an existing RDS via import_genotypes(path)
    # The current implementation of import_genotypes(path) handles "rds" extension.

    obj$import_genotypes(rds_path)

    # Verify it loaded
    expect_equal(nrow(obj$geno), n)

    # Destroy
    rm(obj)
    gc()

    # assert file STILL exists
    expect_true(file.exists(bk_path))
    expect_true(file.exists(rds_path))

    # Cleanup manually
    unlink(bk_path)
    unlink(rds_path)
})
