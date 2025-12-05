## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(HaploGeno)
library(bigstatsr)
library(future)


## ----data_sim-----------------------------------------------------------------
set.seed(123)
n <- 200
m <- 2000

# Generate random genotypes (0, 1, 2)
X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)

# Introduce 1% missing data
na_indices <- sample(length(X), length(X) * 0.01)
X[na_indices] <- NA

# Create a dummy marker map
map <- data.frame(
  chr = rep(1, m),
  id = paste0("snp", 1:m),
  pos = 1:m,
  ref = "A",
  alt = "G"
)

# Simulate a phenotype with effects in the first 100 markers
eff <- rnorm(100)
y <- as.vector(X[, 1:100] %*% eff + rnorm(n))
# Handle NAs in phenotype simulation (simple fix for demo)
y[is.na(y)] <- rnorm(sum(is.na(y)))


## ----init---------------------------------------------------------------------
# Initialize
haplo <- HaploObject$new(tempfile())

# Import genotypes (matrix)
# For PLINK: haplo$import_genotypes("path/to/file.bed")
# For VCF:   haplo$import_genotypes("path/to/file.vcf")
haplo$import_genotypes(X)

# Load map and phenotypes
haplo$load_map(map)
haplo$load_pheno(y)

print(haplo)


## ----preprocessing------------------------------------------------------------
# Filter monomorphic markers
haplo$filter_monomorphic()

# Impute missing genotypes
haplo$impute_genotypes(method = "mean")


## ----blocking-----------------------------------------------------------------
# Define blocks with r2 threshold of 0.1
haplo$define_blocks_ld(r2_threshold = 0.1, window_size = 100)

print(haplo)


## ----encoding-----------------------------------------------------------------
# Use 2 cores for parallel execution
haplo$encode_haplotypes(n_cores = 2)


## ----hrm----------------------------------------------------------------------
haplo$compute_hrm(n_cores = 2)


## ----cv-----------------------------------------------------------------------
# Run 5-fold Cross-Validation
cv_results <- haplo$cross_validate(k = 5, n_cores = 2)

best_lambda <- cv_results$best_lambda
message("Best lambda: ", best_lambda)


## ----fitting------------------------------------------------------------------
alpha <- haplo$fit_krr(lambda = best_lambda)


## ----local_gebv---------------------------------------------------------------
# Estimate marker effects (required for local GEBV)
haplo$estimate_marker_effects(lambda = best_lambda)

# Calculate local GEBVs
haplo$calculate_local_gebv(n_cores = 2)

# Perform significance testing (simple variance test)
haplo$test_significance()


## ----plotting, fig.width=7, fig.height=5--------------------------------------
haplo$plot_manhattan(threshold = 0.05)


## ----saving-------------------------------------------------------------------
haplo$save_project(file.path(tempdir(), "haplo_project.rds"))

