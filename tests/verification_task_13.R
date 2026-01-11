library(testthat)

# Source script
source("standalone_scripts/13_export_hrm_asreml.R")

context("Verification Task: ASReml Export")

test_that("export_hrm_for_asreml generates correct .giv and .map files", {
    # 1. Setup Mock Data
    n_mk <- 10
    n_ind <- 4

    # Genotypes (Haplotype IDs)
    # Ind 1 & 2 Identical -> G=1
    # Ind 3 & 4 Different
    haplo <- matrix(
        c(
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            1, 1, 1, 1, 1, 2, 2, 2, 2, 2
        ),
        nrow = 4, byrow = TRUE
    )
    rownames(haplo) <- paste0("Ind", 1:4)

    # Output locations
    out_prefix <- "test_asreml_export"

    # 2. Run Export
    # Capture output to avoid cluttering console
    capture.output({
        export_hrm_for_asreml(haplo, out_prefix = out_prefix, bending_factor = 0.01)
    })

    # 3. Verify Files Existance
    giv_file <- paste0(out_prefix, ".giv")
    map_file <- paste0(out_prefix, ".map")

    expect_true(file.exists(giv_file))
    expect_true(file.exists(map_file))

    # 4. Verify Content

    # Check Map
    map_df <- read.csv(map_file)
    expect_equal(nrow(map_df), 4)
    expect_equal(map_df$ID, rownames(haplo))
    expect_equal(map_df$Index, 1:4)

    # Check GIV
    # Should have 3 columns: Row, Col, Value
    giv_df <- read.table(giv_file)
    expect_equal(ncol(giv_df), 3)

    # Check Symmetry/Indices
    # Row >= Col usually for lower triangle
    expect_true(all(giv_df[, 1] >= giv_df[, 2]))

    # 5. Cleanup
    file.remove(giv_file)
    file.remove(map_file)
})
