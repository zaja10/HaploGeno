library(testthat)
library(HaploGeno)
library(bigstatsr)

test_that("Safe Imputation rounds values", {
    # Create a small FBM with missing values
    temp_file <- tempfile()
    fbm <- FBM(5, 5, type = "double", backingfile = temp_file)

    # Set some values
    fbm[] <- 0
    # Col 1: 1, 2, NA, 1, 2. Mean = 1.5. Round -> 2.
    fbm[1, 1] <- 1
    fbm[2, 1] <- 2
    fbm[3, 1] <- NA
    fbm[4, 1] <- 1
    fbm[5, 1] <- 2

    # Create HaploObject
    haplo <- HaploObject$new(temp_file)
    haplo$geno <- fbm

    # Impute
    haplo$impute_genotypes(method = "mean")

    # Check value
    val <- fbm[3, 1]
    expect_equal(val, 2) # Should be rounded to 2, not 1.5
})

test_that("Monomorphic markers are filtered", {
    temp_file <- tempfile()
    fbm <- FBM(10, 5, type = "double", backingfile = temp_file)

    # Random data
    fbm[] <- sample(0:2, 50, replace = TRUE)

    # Make col 3 monomorphic (all 0)
    fbm[, 3] <- 0

    haplo <- HaploObject$new(temp_file)
    haplo$geno <- fbm

    haplo$filter_monomorphic()

    expect_equal(length(haplo$active_markers), 4)
    expect_false(3 %in% haplo$active_markers)

    # Check define_haploblocks(method="fixed", uses active markers
    # Mock map
    haplo$map <- data.frame(chr = rep(1, 5), id = paste0("m", 1:5), pos = 1:5, ref = "A", alt = "T")

    haplo$define_haploblocks(method="fixed", window_size = 2)

    # Should have 2 blocks (4 markers / 2 = 2 blocks)
    expect_equal(nrow(haplo$blocks), 2)

    # First block: Start 1, End 2. Indices into active_markers.
    # Active markers: 1, 2, 4, 5.
    # Block 1: 1, 2 -> Real 1, 2.
    # Block 2: 3, 4 -> Real 4, 5.

    # Check encode_haplotypes
    haplo$encode_haplotypes()
    expect_equal(ncol(haplo$haplo_mat), 2)
})

test_that("Dimension consistency checks work", {
    temp_file <- tempfile()
    fbm <- FBM(10, 5, type = "double", backingfile = temp_file)

    haplo <- HaploObject$new(temp_file)
    haplo$geno <- fbm

    # Pheno with 9 samples
    pheno <- rnorm(9)

    expect_error(haplo$load_pheno(pheno), "Number of phenotypes")

    # Pheno with 10 samples
    pheno <- rnorm(10)
    expect_silent(haplo$load_pheno(pheno))
})

test_that("Save and Load Project works", {
    temp_dir <- tempdir()
    proj_name <- "test_proj"
    rds_path <- file.path(temp_dir, paste0(proj_name, ".rds"))
    bk_path <- file.path(temp_dir, proj_name) # FBM backing file

    # Create object
    fbm <- FBM(10, 5, type = "double", backingfile = bk_path)
    fbm[] <- 1

    haplo <- HaploObject$new(bk_path)
    haplo$geno <- fbm
    haplo$pheno <- rnorm(10)

    # Save
    haplo$save_project(rds_path)

    expect_true(file.exists(rds_path))

    # Load
    loaded <- load_haplo_project(rds_path)
    expect_true(!is.null(loaded$geno))
    expect_equal(loaded$geno[1, 1], 1)

    # Test moving (simulated by renaming backing file and trying to load from new location?
    # Or moving both to new dir)

    new_dir <- file.path(temp_dir, "moved")
    dir.create(new_dir)

    file.rename(rds_path, file.path(new_dir, paste0(proj_name, ".rds")))
    file.rename(paste0(bk_path, ".bk"), file.path(new_dir, paste0(proj_name, ".bk")))
    # Note: FBM creates .rds too usually, but we initialized with backingfile path.
    # FBM constructor creates .bk. It also creates .rds if create_bk=TRUE (default).
    # We need to move that too if it exists.
    if (file.exists(paste0(bk_path, ".rds"))) {
        file.rename(paste0(bk_path, ".rds"), file.path(new_dir, paste0(proj_name, ".rds.fbm")))
        # Wait, FBM rds name is usually just .rds. But here we have project .rds.
        # FBM backingfile argument defines the prefix.
    }

    # Load from new location
    new_rds_path <- file.path(new_dir, paste0(proj_name, ".rds"))
    loaded_moved <- load_haplo_project(new_rds_path)

    expect_true(!is.null(loaded_moved$geno))
    # Check if attached
    # Accessing data should work
    expect_equal(loaded_moved$geno[1, 1], 1)
})
