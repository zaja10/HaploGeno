library(testthat)
library(R6)
library(bigstatsr)
library(future)
source("../../R/HaploObject.R")

test_that("Parallel execution works and matches sequential", {
  set.seed(123)
  n <- 50
  m <- 100
  
  geno_mat <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
  map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
  y <- rnorm(n)
  
  tmp_file <- tempfile()
  haplo <- HaploObject$new(tmp_file)
  haplo$import_genotypes(geno_mat)
  haplo$load_map(map)
  haplo$load_pheno(y)
  haplo$define_blocks_fixed(10)
  
  # 1. Parallel Encoding
  # Run sequential first
  haplo$encode_haplotypes(n_cores = 1)
  seq_mat <- haplo$haplo_mat
  
  # Run parallel
  haplo$encode_haplotypes(n_cores = 2)
  par_mat <- haplo$haplo_mat
  
  expect_equal(seq_mat, par_mat)
  
  # 2. Parallel HRM
  haplo$compute_hrm(n_cores = 1)
  seq_hrm <- haplo$hrm
  
  haplo$compute_hrm(n_cores = 2)
  par_hrm <- haplo$hrm
  
  expect_equal(seq_hrm, par_hrm)
  
  # 3. Parallel localGEBV
  haplo$estimate_marker_effects(lambda = 0.1)
  
  haplo$calculate_local_gebv(n_cores = 1)
  seq_gebv <- haplo$local_gebv
  
  haplo$calculate_local_gebv(n_cores = 2)
  par_gebv <- haplo$local_gebv
  
  expect_equal(seq_gebv$matrix, par_gebv$matrix)
  expect_equal(seq_gebv$variances, par_gebv$variances)
  
  message("Parallel execution verified.")
})
