library(testthat)
library(R6)
library(bigstatsr)
source("../../R/HaploObject.R")
source("../../R/S3Methods.R")

test_that("Imputation works", {
  set.seed(123)
  n <- 50
  m <- 50
  
  # Create matrix with NAs
  X <- matrix(rnorm(n * m), n, m)
  X[1, 1] <- NA
  X[10, 10] <- NA
  
  haplo <- HaploObject$new(tempfile())
  haplo$import_genotypes(X)
  
  # Check NAs exist
  expect_true(any(is.na(haplo$geno[])))
  
  # Impute
  haplo$impute_genotypes(method = "mean")
  
  # Check NAs gone
  expect_false(any(is.na(haplo$geno[])))
  
  # Check value (mean of col 1 excluding NA)
  mu <- mean(X[, 1], na.rm = TRUE)
  expect_equal(haplo$geno[1, 1], mu)
})

test_that("Manhattan plot runs", {
  haplo <- HaploObject$new(tempfile())
  
  # Mock data for plotting
  haplo$significance <- data.frame(
    BlockID = 1:10,
    P_Value = runif(10)
  )
  haplo$blocks <- data.table::data.table(
    BlockID = 1:10,
    Start = seq(1, 100, by=10),
    End = seq(10, 100, by=10)
  )
  
  # Should run without error and return data invisibly
  df <- haplo$plot_manhattan()
  expect_true(is.data.frame(df))
  expect_true("logP" %in% names(df))
})

test_that("S3 methods work", {
  haplo <- HaploObject$new(tempfile())
  expect_output(print(haplo), "HaploObject")
  expect_output(summary(haplo), "HaploObject")
})
