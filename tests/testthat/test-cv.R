library(testthat)
library(R6)
library(bigstatsr)
library(future)
source("../../R/HaploObject.R")

test_that("Cross-Validation runs and selects lambda", {
  set.seed(123)
  n <- 50
  m <- 50
  
  X <- matrix(rnorm(n * m), n, m)
  K <- tcrossprod(X) / m
  
  # Simulate y with some signal
  y <- X[, 1] * 2 + rnorm(n)
  
  haplo <- HaploObject$new(tempfile())
  haplo$hrm <- K
  haplo$pheno <- y
  
  # Run CV
  # Use a small grid
  lambdas <- c(10, 1, 0.1, 0.01, 0.001)
  folds <- sample(rep(1:3, length.out = n))
  
  cv_res <- haplo$cross_validate(k = 3, lambdas = lambdas, n_cores = 1, folds = folds)
  
  expect_true(is.list(cv_res))
  expect_true("best_lambda" %in% names(cv_res))
  expect_true("mse" %in% names(cv_res))
  expect_equal(length(cv_res$mse), length(lambdas))
  
  message("Best lambda: ", cv_res$best_lambda)
  
  # Check parallel execution
  cv_res_par <- haplo$cross_validate(k = 3, lambdas = lambdas, n_cores = 2, folds = folds)
  expect_equal(cv_res$best_lambda, cv_res_par$best_lambda)
})
