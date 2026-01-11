library(testthat)
library(ggplot2)
library(data.table)

# Source scripts
source("standalone_scripts/11_vis_haplo_biplot.R")

context("Verification Task: Biplot")

test_that("Script 11 runs and produces a plot", {
    # 1. Setup Mock Data
    n_mk <- 50
    n_ind <- 20

    # Mock GEBVs (Matrix)
    # Create some structure for PCA to find
    # 2 groups, 2 blocks driving diff
    local_gebv <- matrix(rnorm(n_mk * n_ind), nrow = n_ind)
    rownames(local_gebv) <- paste0("Ind", 1:n_ind)

    # Group 1: High on Block 1
    local_gebv[1:10, 1] <- local_gebv[1:10, 1] + 5
    # Group 2: High on Block 2
    local_gebv[11:20, 2] <- local_gebv[11:20, 2] + 5

    # Mock Blocks
    blocks <- data.frame(
        BlockID = 1:n_mk,
        Chr = rep(c("1A", "1B"), each = 25),
        Start = 1:n_mk * 100
    )

    # Mock Groups
    groups <- factor(rep(c("A", "B"), each = 10))

    # 2. Run Plotting (Script 11)
    p <- plot_haplo_biplot(
        local_gebv_list = local_gebv,
        blocks_df = blocks,
        genotype_groups = groups,
        top_vectors = 5
    )

    # Check output
    expect_true(inherits(p, "ggplot"))

    # Save for manual inspection if needed
    # ggsave("test_biplot.png", p)
})
