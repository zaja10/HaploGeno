library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")

test_that("HaploObject data loading works", {
  # Create dummy data
  n <- 10
  m <- 5
  geno <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
  pheno <- rnorm(n)

  # Initialize object
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)

  # Test import_genotypes (matrix)
  haplo$import_genotypes(geno)
  expect_true(!is.null(haplo$geno))
  expect_equal(nrow(haplo$geno), n)
  expect_equal(ncol(haplo$geno), m)

  # Test load_map
  haplo$load_map(map)
  expect_true(!is.null(haplo$map))
  expect_equal(nrow(haplo$map), m)

  # Test load_pheno
  haplo$load_pheno(pheno)
  expect_true(!is.null(haplo$pheno))
  expect_equal(length(haplo$pheno), n)

  # Test validations
  expect_warning(haplo$load_map(rbind(map, map)), "Number of markers")
  expect_error(haplo$load_pheno(c(pheno, 1)), "Number of phenotypes")
})
