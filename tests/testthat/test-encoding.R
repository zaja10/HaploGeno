library(testthat)
library(R6)
library(bigstatsr)
library(Rcpp)

# Source R file


test_that("HaploObject blocking and encoding works", {
  # Create dummy data
  n <- 10
  m <- 20
  geno <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")

  # Initialize object
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)
  haplo$import_genotypes(geno)
  haplo$load_map(map)

  # Test define_haploblocks(method="fixed",
  window_size <- 5
  haplo$define_haploblocks(method = "fixed", window_size = window_size)

  expect_true(!is.null(haplo$blocks))
  expect_equal(nrow(haplo$blocks), m / window_size)
  expect_equal(haplo$blocks$Start[1], 1)
  expect_equal(haplo$blocks$End[1], 5)

  # Test encode_haplotypes
  haplo$encode_haplotypes()

  expect_true(!is.null(haplo$haplo_mat))
  expect_equal(nrow(haplo$haplo_mat), n)
  expect_equal(ncol(haplo$haplo_mat), nrow(haplo$blocks))

  # Check uniqueness logic manually for first block
  block1 <- geno[, 1:5]
  # Create string representation
  strs <- apply(block1, 1, paste, collapse = "")
  unique_strs <- unique(strs)
  # Number of unique haplotypes should match max ID in that column (if IDs are sequential)
  # Note: IDs might not be 1..K if not re-indexed, but our C++ does 1..K
  expect_equal(length(unique(haplo$haplo_mat[, 1])), length(unique_strs))
})
