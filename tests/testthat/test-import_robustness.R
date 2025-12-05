test_that("import_genotypes handles dosage CSV", {
    # Create temporary CSV
    tf <- tempfile(fileext = ".csv")

    # 5 samples, 5 markers
    # First col = ID
    df <- data.frame(
        ID = paste0("Ind", 1:5),
        M1 = c(0, 1, 2, 0, 1),
        M2 = c(2, 2, 0, 1, 0),
        M3 = c(0, 0, 0, 0, 0),
        M4 = c(1, 1, 1, 1, 1),
        M5 = c(2, 1, 0, 1, 2)
    )

    write.csv(df, tf, row.names = FALSE)

    # Init object
    ho <- HaploObject$new(tempfile())

    # Import
    expect_error(ho$import_genotypes(tf), NA)

    # Check dimensions
    expect_equal(nrow(ho$geno), 5)
    expect_equal(ncol(ho$geno), 5)

    # Check values (M1)
    # FBM is column-major usually, but access via [] works
    # We removed ID column, so M1 is first column
    expect_equal(ho$geno[, 1], df$M1)

    # Check map
    expect_true(!is.null(ho$map))
    expect_equal(nrow(ho$map), 5)
    expect_equal(ho$map$marker, paste0("M", 1:5))
})

test_that("import_genotypes handles phased/slashed CSV", {
    tf <- tempfile(fileext = ".csv")

    # Phased data
    df <- data.frame(
        ID = paste0("Ind", 1:3),
        M1 = c("0/0", "0/1", "1/1"),
        M2 = c("0|0", "0|1", "1|1"),
        M3 = c("A/A", "A/B", "B/B")
    )

    write.csv(df, tf, row.names = FALSE)

    ho <- HaploObject$new(tempfile())
    ho$import_genotypes(tf)

    expect_equal(nrow(ho$geno), 3)
    expect_equal(ncol(ho$geno), 3)

    # Check conversions
    # 0/0 -> 0
    # 0/1 -> 1
    # 1/1 -> 2
    expect_equal(ho$geno[, 1], c(0, 1, 2))

    # 0|0 -> 0
    # 0|1 -> 1
    # 1|1 -> 2
    expect_equal(ho$geno[, 2], c(0, 1, 2))

    # A/A -> 0
    # A/B -> 1
    # B/B -> 2
    expect_equal(ho$geno[, 3], c(0, 1, 2))
})

test_that("import_genotypes handles missing values", {
    tf <- tempfile(fileext = ".csv")

    # Missing data
    df <- data.frame(
        ID = paste0("Ind", 1:3),
        M1 = c("0", "1", NA),
        M2 = c("0/0", ".", "1/1")
    )

    # Write with NA as empty string or NA
    write.csv(df, tf, row.names = FALSE, na = ".")

    ho <- HaploObject$new(tempfile())

    # Expect warning about missing values
    expect_warning(ho$import_genotypes(tf), "missing values")

    expect_true(is.na(ho$geno[3, 1]))
    expect_true(is.na(ho$geno[2, 2]))
})
