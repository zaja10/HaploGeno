# Script: 09_vis_association.R
# Description: Tools for interpreting association results. 
#              Includes Volcano Plots and Haplotype Stacking Analysis (Identification & Scoring).
# Dependencies: graphics, stats, data.table

library(graphics)
library(stats)
library(data.table)

# --- 1. Stacking Analysis Functions ---

#' Identify Superior Haplotypes
#'
#' Finds the specific haplotype allele with the highest positive effect 
#' in the most significant blocks.
#'
#' @param local_gebv_list Output from Script 4 (must contain $matrix).
#' @param haplo_mat Output from Script 2 (Haplotype IDs).
#' @param significance_df Output from Script 5 (P-values).
#' @param top_n Number of top significant blocks to use (default: 50).
#' @return A data.frame defining the best haplotype for each top block.
identify_superior_haplotypes_func <- function(local_gebv_list, haplo_mat, significance_df, top_n = 50) {
    
    # Input validation / Unpacking
    if (is.list(local_gebv_list) && "matrix" %in% names(local_gebv_list)) {
        gebv_mat <- local_gebv_list$matrix
    } else {
        gebv_mat <- as.matrix(local_gebv_list)
    }
    
    # Sort significance table by P-value (ascending)
    # Ensure we only consider valid blocks that exist in our matrices
    # (Assuming BlockIDs align with column indices 1..M for simplicity, or match names)
    
    # Robust sort
    sorted_sig <- significance_df[order(significance_df$P_Value), ]
    
    # Select top N blocks
    top_blocks <- head(sorted_sig, top_n)
    
    message("Identifying superior haplotypes for top ", nrow(top_blocks), " blocks...")
    
    results <- list()
    
    for (i in 1:nrow(top_blocks)) {
        block_id <- top_blocks$BlockID[i]
        
        # Get data for this block
        # Assuming haplo_mat and gebv_mat columns correspond to BlockIDs 1..M
        # If columns are named "B1", "B2", we rely on index if BlockID is integer index
        
        # Check index bounds
        if (block_id > ncol(gebv_mat) || block_id > ncol(haplo_mat)) {
            warning("BlockID ", block_id, " exceeds matrix dimensions. Skipping.")
            next
        }
        
        gebvs <- gebv_mat[, block_id]
        haplos <- haplo_mat[, block_id]
        
        # Calculate mean GEBV for each haplotype allele
        # aggregate returns data.frame(Group, x)
        agg <- aggregate(gebvs, by = list(haplos), FUN = mean)
        colnames(agg) <- c("HaploID", "MeanGEBV")
        
        # Identify the "Best" haplotype (Highest Positive Effect)
        best_idx <- which.max(agg$MeanGEBV)
        best_haplo <- agg$HaploID[best_idx]
        effect_size <- agg$MeanGEBV[best_idx]
        
        results[[i]] <- data.frame(
            BlockID = block_id,
            BestHaploID = best_haplo,
            EffectSize = effect_size
        )
    }
    
    return(do.call(rbind, results))
}

#' Calculate Stacking Scores
#'
#' Scores individuals based on how many "Superior Haplotypes" they possess.
#'
#' @param haplo_mat Matrix of Haplotype IDs.
#' @param superior_df Data frame returned by identify_superior_haplotypes_func.
#' @return Integer vector of scores (0 to top_n).
score_stacking_func <- function(haplo_mat, superior_df) {
    message("Scoring individuals based on ", nrow(superior_df), " superior haplotypes...")
    
    n_samples <- nrow(haplo_mat)
    scores <- integer(n_samples) # Initialize with 0
    
    for (i in 1:nrow(superior_df)) {
        b_id <- superior_df$BlockID[i]
        target_allele <- superior_df$BestHaploID[i]
        
        # Check which individuals have the target allele at this block
        # Returns logical vector, convert to integer (0/1)
        matches <- (haplo_mat[, b_id] == target_allele)
        
        scores <- scores + as.integer(matches)
    }
    
    return(scores)
}

# --- 2. Visualization Functions ---

#' Volcano Plot of Block Significance
#'
#' @param significance_df Data frame with 'Variance' and 'P_Value' columns.
#' @param p_threshold P-value threshold to highlight points (default 0.05).
#' @param main Title of the plot.
plot_volcano <- function(significance_df, p_threshold = 0.05, main = "HaploBlock Volcano Plot") {
    if (missing(significance_df) || is.null(significance_df)) stop("significance_df is required.")

    df <- significance_df

    # Safety: Ensure no log(0)
    safe_p <- df$P_Value
    safe_p[safe_p <= 0] <- 1e-300 # Tiny float
    
    y <- -log10(safe_p)
    x <- df$Variance

    # Highlight Logic
    is_sig <- df$P_Value < p_threshold
    
    # Base Plot (Transparent Grey for non-significant)
    # Type 'n' initializes plot without points
    plot(x, y, type = "n",
        xlab = "Local Genetic Variance",
        ylab = expression(-log[10](italic(p))),
        main = main
    )
    
    grid()
    
    # Plot non-significant points
    points(x[!is_sig], y[!is_sig], pch = 19, col = rgb(0.5, 0.5, 0.5, 0.3))
    
    # Plot significant points (Red)
    if (any(is_sig)) {
        points(x[is_sig], y[is_sig], pch = 19, col = rgb(1, 0, 0, 0.6), cex = 1.0)
    }
    
    # Threshold Line
    abline(h = -log10(p_threshold), col = "blue", lty = 2, lwd = 1.5)
    text(x = min(x), y = -log10(p_threshold), labels = paste0("p < ", p_threshold), 
         pos = 3, col = "blue", cex = 0.8)
}

#' Plot Stacking Trend
#'
#' Visualizes the relationship between the Stacking Score and Phenotype.
#'
#' @param scores Integer vector of stacking scores.
#' @param pheno_vec Numeric vector of phenotypes.
#' @param main Title of the plot.
plot_stacking_trend <- function(scores, pheno_vec, main = "Stacking Validation") {
    if (missing(scores) || is.null(scores)) stop("scores are required.")
    if (missing(pheno_vec) || is.null(pheno_vec)) stop("pheno_vec is required.")
    
    # Align vectors (if named)
    if (!is.null(names(scores)) && !is.null(names(pheno_vec))) {
        common <- intersect(names(scores), names(pheno_vec))
        scores <- scores[common]
        pheno_vec <- pheno_vec[common]
    }
    
    if (length(scores) != length(pheno_vec)) stop("Length of scores and phenotypes must match.")

    # Plot
    plot(scores, pheno_vec,
        xlab = "Stacking Score (N Superior Haplotypes)",
        ylab = "Phenotype",
        main = main,
        pch = 19,
        col = rgb(0, 0, 0, 0.4), # Transparency helps see density
        las = 1
    )

    # Regression line
    fit <- lm(pheno_vec ~ scores)
    abline(fit, col = "red", lwd = 2)

    # Stats
    summ <- summary(fit)
    r2 <- summ$r.squared
    pval <- summ$coefficients[2, 4] # P-value of slope
    
    legend("topleft", 
           legend = c(
               paste0("R2 = ", round(r2, 3)),
               paste0("P = ", format(pval, digits = 3))
           ), 
           bty = "n", text.col = "red")
}
