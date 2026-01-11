library(testthat)
library(ggplot2)
library(data.table)

# Source scripts
source("standalone_scripts/01_haplotype_blocking.R")
source("standalone_scripts/10_vis_chromosome_painting.R")

context("Verification Task: Chromosome Painting")

test_that("Script 10 runs and produces a plot", {
    # 1. Setup Mock Data
    n_mk <- 100
    n_ind <- 5

    geno <- matrix(sample(0:2, n_mk * n_ind, replace = TRUE), nrow = n_ind)
    rownames(geno) <- paste0("Ind", 1:n_ind)

    map <- data.frame(
        chr = rep(c("1A", "1B"), each = 50),
        pos = seq(1, 1000, length.out = 100)
    )

    # 2. Run Block Definition (Script 01)
    # This should now include 'Chr' in blocks_df
    blocks <- define_blocks_ld(geno, map, verbose = FALSE)

    expect_true("Chr" %in% names(blocks))
    expect_true("1A" %in% blocks$Chr)

    # 3. Setup Mock GEBVs (Script 04 output style)
    # num blocks
    n_blocks <- nrow(blocks)
    local_gebv <- matrix(rnorm(n_blocks * n_ind), nrow = n_ind)
    rownames(local_gebv) <- rownames(geno)

    # 4. Setup Mock Genes
    genes <- data.frame(
        Gene = c("GeneA", "GeneB"),
        Chr = c("1A", "1B"),
        Pos = c(250, 750),
        Trait_Category = "Yield"
    )

    # 5. Run Plotting (Script 10)
    p <- plot_chromosome_painting(
        local_gebv_list = local_gebv,
        blocks_df = blocks,
        target_chr = "1A",
        target_genotypes = c("Ind1", "Ind2"),
        highlight_genes_df = genes,
        map_data = map
    )

    # Check output
    expect_true(inherits(p, "ggplot"))

    # Save for manual inspection if needed (optional)
    # ggsave("test_painting.png", p)
})
