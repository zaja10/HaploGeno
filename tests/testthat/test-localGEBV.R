library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")

test_that("localGEBV workflow runs correctly", {
  set.seed(123)
  n <- 100
  m <- 100
  
  # Simulate genotypes with some LD structure
  # Create blocks of 10 markers with high correlation
  geno_mat <- matrix(0, n, m)
  for (i in seq(1, m, by = 10)) {
    # Base vector
    base <- sample(0:2, n, replace = TRUE)
    for (j in 0:9) {
      if (i + j <= m) {
        # High correlation with base
        noise <- sample(c(0, 0, 1), n, replace = TRUE) # Sparse noise
        g <- base + noise
        g[g > 2] <- 2
        geno_mat[, i + j] <- g
      }
    }
  }
  
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
  
  # Simulate phenotype: QTL in block 2 (markers 11-20)
  # Effect is on marker 15
  y <- geno_mat[, 15] * 50 + rnorm(n)
  
  # Initialize
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)
  haplo$import_genotypes(geno_mat)
  haplo$load_map(map)
  haplo$load_pheno(y)
  
  # Step 1: LD Blocking
  # Expect roughly m/10 blocks if LD is perfect, but noise might break them
  haplo$define_haploblocks(method="ld", r2_threshold = 0.5, tolerance = 2)
  expect_true(!is.null(haplo$blocks))
  expect_gt(nrow(haplo$blocks), 1)
  
  # Step 2: Marker Effects
  haplo$estimate_marker_effects(lambda = 0.1)
  expect_true(!is.null(haplo$marker_effects))
  expect_equal(length(haplo$marker_effects), m)
  
  # Step 3 & 4: localGEBV
  haplo$calculate_local_gebv()
  expect_true(!is.null(haplo$local_gebv))
  expect_equal(dim(haplo$local_gebv$matrix), c(n, nrow(haplo$blocks)))
  
  # Step 5: Significance
  sig <- haplo$test_significance()
  expect_true(!is.null(sig))
  expect_true("P_Value" %in% names(sig))
  
  # Check if the block containing marker 15 is significant
  # Find block containing marker 15
  block_idx <- which(haplo$blocks$Start <= 15 & haplo$blocks$End >= 15)
  expect_true(length(block_idx) > 0)
  
  # The p-value for this block should be small (significant)
  p_val <- sig$P_Value[block_idx]
  var_val <- sig$Variance[block_idx]
  
  message("Significant block Variance: ", var_val)
  message("Significant block P-value: ", p_val)
  
  expect_lt(p_val, 0.2)
})
