test_that("Analysis features work correctly", {
  # Create a small mock HaploObject
  # We need a small FBM

  # 50 individuals, 100 markers
  n <- 50
  m <- 100

  # Create temporary backing file
  tmp_file <- tempfile()

  # Random genotypes (0, 1, 2)
  X <- matrix(sample(0:2, n * m, replace = TRUE), n, m)

  # Initialize object
  obj <- HaploObject$new(tmp_file)
  obj$import_genotypes(X)

  # Create mock map
  map <- data.frame(
    chr = rep(1:2, each = 50),
    id = paste0("rs", 1:m),
    pos = 1:m,
    ref = "A",
    alt = "T"
  )
  obj$load_map(map)

  # Create mock phenotype
  # y = X[, 1] * 2 + noise
  y <- X[, 1] * 2 + rnorm(n)
  obj$load_pheno(y)

  # Run pipeline steps needed for analysis
  obj$define_haploblocks(method = "fixed", window_size = 10)
  obj$encode_haplotypes()
  obj$compute_hrm()
  obj$estimate_marker_effects()
  obj$calculate_local_gebv()
  obj$test_significance()

  # Test plot_manhattan
  # We can't easily test the plot output, but we can check if it runs without error
  # and returns the invisible data frame
  p <- obj$plot_manhattan()
  expect_true(inherits(p, "ggplot"))
  df <- p$data
  # Check for essential columns in the plot data
  expect_true(all(c("BlockID", "P_Value", "Pos") %in% names(df)))
  expect_equal(nrow(df), nrow(obj$blocks))

  # Test plot_pca
  # Should run without error
  expect_error(obj$plot_pca(), NA)

  # Test plot_gebv_image
  expect_error(obj$plot_gebv_image(), NA)

  # Test identify_superior_haplotypes
  top_n <- 5
  res <- obj$identify_superior_haplotypes(top_n = top_n)
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), top_n)
  expect_true(all(c("BlockID", "BestHaploID", "EffectSize") %in% names(res)))

  # Test score_stacking
  scores <- obj$score_stacking(res)
  expect_true(is.numeric(scores))
  expect_equal(length(scores), n)
  expect_true(all(scores >= 0))
  expect_true(all(scores <= top_n))

  # Test plot_stacking_trend
  # Should run without error
  expect_error(obj$plot_stacking_trend(scores = scores), NA)
  expect_error(obj$plot_stacking_trend(superior_haplos = res), NA)

  # Clean up
  unlink(paste0(tmp_file, ".bk"))
  unlink(paste0(tmp_file, ".rds"))
})
