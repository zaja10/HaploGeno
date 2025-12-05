## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(HaploGeno)


## ----load_demo----------------------------------------------------------------
# Load the pre-processed demo dataset
haplo <- load_demo_data()

# Print the object summary to verify data status
print(haplo)


## ----pca, fig.width=6, fig.height=6-------------------------------------------
# Plot the first two principal components of the HRM
haplo$plot_pca()


## ----manhattan, fig.width=7, fig.height=5-------------------------------------
# Plot significance of local GEBV variances
# The threshold defaults to p < 0.05
haplo$plot_manhattan(threshold = 0.05)


## ----heatmap, fig.width=7, fig.height=6---------------------------------------
# Visualize the Local GEBV matrix
# X-axis: Block Index
# Y-axis: Individual Index
# Color: Blue (Low Value) to Red (High Value)
haplo$plot_gebv_image()


## ----superior-----------------------------------------------------------------
# Identify the best haplotype for the top 10 most significant blocks
superior_haplos <- haplo$identify_superior_haplotypes(top_n = 10)

# View the results
print(head(superior_haplos))


## ----scoring------------------------------------------------------------------
# Calculate stacking scores for all individuals
scores <- haplo$score_stacking(superior_haplos)

# View scores for the first few individuals
head(scores)


## ----trend, fig.width=6, fig.height=6-----------------------------------------
# Plot the stacking trend with a regression line and R-squared value
haplo$plot_stacking_trend(scores = scores)

