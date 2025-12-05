test_that("plot_haplo_biplot works", {
    # Initialize object
    tmp_file <- tempfile()
    haplo <- HaploObject$new(tmp_file)

    # Mock necessary data
    n_ind <- 20
    n_blocks <- 10

    # HRM (Identity for simplicity)
    haplo$hrm <- diag(n_ind)

    # Local GEBVs (Random)
    local_gebv <- matrix(rnorm(n_ind * n_blocks), nrow = n_ind, ncol = n_blocks)
    haplo$local_gebv <- list(matrix = local_gebv)

    # Significance table
    haplo$significance <- data.table::data.table(
        BlockID = 1:n_blocks,
        Variance = runif(n_blocks),
        P_Value = runif(n_blocks)
    )

    # Test basic run
    expect_error(res <- haplo$plot_haplo_biplot(top_n = 5, label_blocks = FALSE), NA)

    # Check result structure
    expect_type(res, "list")
    expect_named(res, c("scores", "loadings"))
    expect_equal(nrow(res$loadings), 5)
    expect_equal(nrow(res$scores), n_ind)

    # Test validation
    haplo$hrm <- NULL
    expect_error(haplo$plot_haplo_biplot(), "HRM not computed")
})
