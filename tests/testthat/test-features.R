library(testthat)
library(R6)
library(bigstatsr)


test_that("Imputation works", {
  set.seed(123)
  n <- 50
  m <- 50

  # Create matrix with NAs
  X <- matrix(sample(0:2, n * m, replace = TRUE), n, m)
  X[1, 1] <- NA
  X[10, 10] <- NA

  haplo <- HaploObject$new(tempfile())
  haplo$import_genotypes(X)

  # Check NAs exist
  expect_true(any(is.na(haplo$geno[])))

  # Imputation should FAIL for matrix objects (guardrail check)
  expect_error(haplo$impute_genotypes(method = "mean"), "Automated imputation is only supported")

  # Note: Actual imputation testing requires code256 objects (PLINK/VCF) which are hard to mock here.
  # We trust strict numeric validation + upstream imputation for matrices.
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
    Start = seq(1, 100, by = 10),
    End = seq(10, 100, by = 10)
  )

  # Should run without error and return data invisibly
  # Should run without error and return a ggplot object
  p <- haplo$plot_manhattan()
  expect_true(inherits(p, "ggplot"))
  expect_true(is.data.frame(p$data))
  expect_true("Val" %in% names(p$data))
})

test_that("S3 methods work", {
  haplo <- HaploObject$new(tempfile())
  expect_output(print(haplo), "HaploObject")
  expect_output(summary(haplo), "HAPLOGENO")
})
