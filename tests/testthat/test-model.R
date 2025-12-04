library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")

test_that("HaploObject modeling works", {
  # Create dummy data
  n <- 50
  m <- 20
  geno <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
  
  # Simulate phenotype: y = X * beta + e
  # Let's say first 2 blocks have effect
  beta_true <- c(1, -1, rep(0, m/5 - 2)) # Assuming window_size=5 -> 4 blocks
  
  # Initialize object
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)
  haplo$import_genotypes(geno)
  haplo$load_map(map)
  
  # Define blocks and encode
  haplo$define_blocks_fixed(5)
  haplo$encode_haplotypes()
  
  # Simulate phenotype based on haplotypes
  # We need to know which haplotype ID corresponds to what
  # For simplicity, just random phenotype
  pheno <- rnorm(n)
  haplo$load_pheno(pheno)
  
  # Compute HRM
  haplo$compute_hrm()
  expect_true(!is.null(haplo$hrm))
  expect_equal(dim(haplo$hrm), c(n, n))
  
  # Fit KRR
  alpha_hat <- haplo$fit_krr(lambda = 0.1)
  
  expect_true(!is.null(alpha_hat))
  expect_equal(length(alpha_hat), n) # Dual coefficients, one per sample
  
  # Check if it runs without error
})
