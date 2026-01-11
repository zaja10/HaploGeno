# Script: 07_vis_genotype_structure.R
# Description: Visualization functions for Genotype and Haplotype Structure.
#              Includes computation of the Haplotype Relationship Matrix (HRM).
# Dependencies: graphics, stats, future.apply (for HRM computation)

library(graphics)
library(stats)
library(future.apply)

# --- 1. Computation Helper ---

#' Compute Haplotype Relationship Matrix (HRM)
#'
#' Calculates the proportion of shared haplotype blocks between individuals.
#' HRM_ij = (Number of blocks where Ind_i and Ind_j have same allele) / Total Blocks
#'
#' @param haplo_mat Matrix of Haplotype IDs (N Individuals x M Blocks).
#' @param n_cores Number of cores for parallel processing.
#' @return A symmetric N x N matrix (0 to 1).
compute_hrm_func <- function(haplo_mat, n_cores = 1) {
    n_samples <- nrow(haplo_mat)
    n_blocks <- ncol(haplo_mat)
    
    message("Computing HRM for ", n_samples, " individuals across ", n_blocks, " blocks...")
    
    # Chunking for parallel processing
    if (n_cores > 1) {
        future::plan(future::multisession, workers = n_cores)
        chunk_indices <- split(1:n_blocks, cut(1:n_blocks, n_cores, labels = FALSE))
    } else {
        future::plan(future::sequential)
        chunk_indices <- list(1:n_blocks)
    }

    # Map: Compute partial K for each chunk of blocks
    partial_Ks <- future_lapply(chunk_indices, function(indices) {
        K_partial <- matrix(0, n_samples, n_samples)
        
        # Iterate through assigned blocks
        for (b in indices) {
            ids <- haplo_mat[, b]
            # Outer product equality: 1 if match, 0 if not
            # This is O(N^2) per block
            K_partial <- K_partial + outer(ids, ids, `==`)
        }
        return(K_partial)
    }, future.seed = TRUE)

    # Reset plan
    future::plan(future::sequential)

    # Reduce: Sum partial matrices
    K_total <- Reduce(`+`, partial_Ks)

    # Normalize by total number of blocks
    HRM <- K_total / n_blocks
    
    # Add names if available
    if (!is.null(rownames(haplo_mat))) {
        rownames(HRM) <- colnames(HRM) <- rownames(haplo_mat)
    }

    return(HRM)
}

# --- 2. Visualization Functions ---

#' Plot Distribution of Haplotype Block Sizes
#'
#' @param blocks_df Data frame with columns Start and End.
#' @param main Title of the plot.
#' @param col Color of histogram bars.
plot_block_sizes <- function(blocks_df, main = "Distribution of Haplotype Block Sizes", col = "steelblue") {
    if (missing(blocks_df) || is.null(blocks_df)) stop("blocks_df is required.")

    # Calculate size in markers
    sizes <- blocks_df$End - blocks_df$Start + 1

    hist(sizes,
        breaks = 30,
        col = col,
        border = "white",
        main = main,
        xlab = "Number of Markers per Block"
    )

    # Add median line
    med_val <- median(sizes)
    abline(v = med_val, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median:", med_val), col = "red", lty = 2, bty = "n")
}

#' Plot Haplotype Relationship Matrix (HRM)
#'
#' @param hrm_mat Square numeric matrix (HRM).
#' @param main Title of the plot.
plot_hrm <- function(hrm_mat, main = "Haplotype Relationship Matrix") {
    if (missing(hrm_mat) || is.null(hrm_mat)) stop("hrm_mat is required. Run compute_hrm_func() first.")

    # Heatmap palette (White -> Orange -> Dark Red)
    pal <- colorRampPalette(c("white", "orange", "darkred"))(100)

    # Image plot nuances:
    # 1. image() puts (0,0) at bottom-left. 
    # 2. To visualize matrix 'properly' (row 1 at top), we flip the Y-axis data.
    
    image(1:nrow(hrm_mat), 1:ncol(hrm_mat), t(hrm_mat)[, nrow(hrm_mat):1],
        col = pal,
        axes = FALSE,
        main = main,
        xlab = "Individuals",
        ylab = "Individuals"
    )

    # Add box around plot
    box()
}

#' Plot PCA of Haplotype Relationship Matrix
#'
#' @param hrm_mat Square numeric matrix (HRM).
#' @param groups Optional vector of groups for coloring (length = nrow(hrm_mat)).
#' @param main Title of the plot.
plot_pca <- function(hrm_mat, groups = NULL, main = "PCA of Haplotype Relationship Matrix") {
    if (missing(hrm_mat) || is.null(hrm_mat)) stop("hrm_mat is required.")
    
    # Safety Check: Matrix must be symmetric
    if (!isSymmetric(hrm_mat)) {
        warning("Matrix is not symmetric. Using symmetric part for Eigen decomposition.")
        hrm_mat <- (hrm_mat + t(hrm_mat)) / 2
    }

    message("Calculating Eigen decomposition of HRM...")
    eig <- eigen(hrm_mat, symmetric = TRUE)

    # Extract top 2 PCs
    pc1 <- eig$vectors[, 1]
    pc2 <- eig$vectors[, 2]

    # Variance explained (Eigenvalues / Sum of Eigenvalues)
    var_expl <- eig$values / sum(eig$values) * 100
    
    # Handle small negative eigenvalues (numerical noise)
    var_expl <- pmax(var_expl, 0)

    # Plot Labels
    xlab <- paste0("PC1 (", round(var_expl[1], 1), "%)")
    ylab <- paste0("PC2 (", round(var_expl[2], 1), "%)")

    # Group Coloring Logic
    col_vec <- "black"
    pch_vec <- 19
    legend_text <- NULL
    
    if (!is.null(groups)) {
        if (length(groups) != nrow(hrm_mat)) {
            warning("Length of groups does not match number of individuals. Ignoring groups.")
        } else {
            groups <- as.factor(groups)
            # Use a slightly nicer palette than integer defaults if possible
            if (length(levels(groups)) <= 8) {
                cols <- c("red", "blue", "green", "purple", "orange", "brown", "cyan", "magenta")
                col_vec <- cols[as.integer(groups)]
            } else {
                col_vec <- as.integer(groups)
            }
            legend_text <- levels(groups)
        }
    }

    plot(pc1, pc2,
        xlab = xlab, ylab = ylab,
        main = main,
        pch = pch_vec,
        col = col_vec,
        bty = "L" # L-shaped box
    )
    
    grid()

    if (!is.null(legend_text)) {
        legend("topright", legend = legend_text, col = unique(col_vec), pch = 19, bty = "n")
    }
}
