# Script: 09_vis_association.R
# Description: Visualization functions for Association Results (Volcano Plot, Stacking).
# Dependencies: graphics (base R)

#' Volcano Plot of Block Significance
#'
#' @param significance_df Data frame with 'Variance' and 'P_Value' columns.
#' @param p_threshold P-value threshold to highlight points (default 0.05).
#' @param main Title of the plot.
plot_volcano <- function(significance_df, p_threshold = 0.05, main = "HaploBlock Volcano Plot") {
    if (missing(significance_df) || is.null(significance_df)) stop("significance_df is required.")

    df <- significance_df

    y <- -log10(df$P_Value)
    x <- df$Variance

    # Highlight color
    # Note: y is -log10(p), so check if p < threshold => y > -log10(threshold)
    col_vec <- ifelse(df$P_Value < p_threshold, "red", "black")

    plot(x, y,
        pch = 19,
        col = rgb(0, 0, 0, 0.3), # Base transparency for all
        xlab = "Local Genetic Variance",
        ylab = "-log10(P-Value)",
        main = main
    )

    # Overlay significant hits with solid color
    sig_idx <- which(df$P_Value < p_threshold)
    if (length(sig_idx) > 0) {
        points(x[sig_idx], y[sig_idx], col = "red", pch = 19, cex = 0.8)
    }

    grid()
}

#' Plot Stacking Trend
#'
#' Visualizes the relationship between the Stacking Score (count of superior haplotypes) and Phenotype.
#'
#' @param scores Integer vector of stacking scores.
#' @param pheno_vec Numeric vector of phenotypes.
#' @param main Title of the plot.
plot_stacking_trend <- function(scores, pheno_vec, main = "Stacking Validation") {
    if (missing(scores) || is.null(scores)) stop("scores are required.")
    if (missing(pheno_vec) || is.null(pheno_vec)) stop("pheno_vec is required.")
    if (length(scores) != length(pheno_vec)) stop("Length of scores and phenotypes must match.")

    plot(scores, pheno_vec,
        xlab = "Stacking Score (N Superior Haplotypes)",
        ylab = "Phenotype",
        main = main,
        pch = 19,
        col = rgb(0, 0, 0, 0.5)
    )

    # Regression line
    fit <- lm(pheno_vec ~ scores)
    abline(fit, col = "red", lwd = 2)

    # Add R2
    r2 <- summary(fit)$r.squared
    legend("topleft", legend = paste0("R2 = ", round(r2, 3)), bty = "n")
}
