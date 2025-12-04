library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")

test_that("KRR captures haplotype effects in simulation", {
  set.seed(42)
  n <- 100
  m <- 200
  
  # Simulate genotypes
  geno <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
  
  # Initialize object
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)
  haplo$import_genotypes(geno)
  haplo$load_map(map)
  
  # Define blocks (size 10)
  haplo$define_blocks_fixed(10)
  haplo$encode_haplotypes()
  
  # Simulate Phenotype based on Haplotypes
  # Let's say the first block has a strong effect
  # We need to access the haplotype IDs to simulate the phenotype
  haplo_ids <- haplo$haplo_mat[, 1]
  
  # Assign effects to specific haplotype IDs in block 1
  # e.g. ID 1 has effect +2, ID 2 has effect -2
  true_effects <- rnorm(max(haplo_ids))
  y_true <- true_effects[haplo_ids]
  
  # Add some noise
  y <- y_true + rnorm(n, sd = 0.5)
  
  haplo$load_pheno(y)
  
  # Fit KRR
  haplo$compute_hrm()
  alpha <- haplo$fit_krr(lambda = 0.1)
  
  # Predict on training data: y_hat = K * alpha
  y_hat <- haplo$hrm %*% alpha
  
  # Check correlation
  cor_val <- cor(y, y_hat)
  
  expect_gt(cor_val, 0.5)
  message("Prediction Correlation: ", round(cor_val, 3))
})
