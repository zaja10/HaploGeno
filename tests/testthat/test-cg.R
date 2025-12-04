library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")

test_that("CG solver matches direct solver", {
  set.seed(123)
  n <- 100
  m <- 50
  
  # Create a random positive definite matrix K
  X <- matrix(rnorm(n * m), n, m)
  K <- tcrossprod(X)
  K <- K / m
  
  y <- rnorm(n)
  
  # Create mock HaploObject
  haplo <- HaploObject$new(tempfile())
  haplo$hrm <- K
  haplo$pheno <- y
  
  # 1. Direct Solver
  alpha_direct <- haplo$fit_krr(lambda = 0.1, use_cg = FALSE)
  
  # 2. CG Solver
  alpha_cg <- haplo$fit_krr(lambda = 0.1, use_cg = TRUE, tol = 1e-6)
  
  # Compare
  diff <- max(abs(alpha_direct - alpha_cg))
  message("Max difference between Direct and CG: ", diff)
  
  expect_lt(diff, 1e-4)
})
