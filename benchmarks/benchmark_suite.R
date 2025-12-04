library(bigstatsr)
library(future)
library(rrBLUP)
library(BGLR)
library(microbenchmark)

# Source package files directly to avoid build issues
source("R/HaploObject.R")
source("R/S3Methods.R")

# Benchmark function
run_benchmark <- function(n, m) {
    message("Benchmarking with N=", n, ", M=", m, "...")
    
    # Simulate data
    set.seed(123)
    X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
    map <- data.frame(chr = rep(1, m), id = paste0("snp", 1:m), pos = 1:m, ref = "A", alt = "G")
    eff <- rnorm(m / 10)
    y <- as.vector(X[, 1:(m/10)] %*% eff + rnorm(n))
    
    # rrBLUP setup
    # rrBLUP expects markers in (-1, 0, 1)
    X_rr <- X - 1
    
    # HaploGeno setup
    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(X)
    haplo$load_map(map)
    haplo$load_pheno(y)
    
    res <- microbenchmark(
        HaploGeno = {
            haplo$define_blocks_ld(r2_threshold = 0.1, window_size = 100)
            haplo$encode_haplotypes(n_cores = 1) 
            haplo$compute_hrm(n_cores = 1)
            haplo$fit_krr(lambda = 1.0)
        },
        rrBLUP = {
            K <- A.mat(X_rr)
            mixed.solve(y, K = K)
        },
        times = 5
    )
    
    return(res)
}

# Run benchmarks
results_500 <- run_benchmark(500, 2000)
print(results_500)

# Plot using base R
boxplot(results_500, main = "HaploGeno vs rrBLUP (N=500, M=2000)")
