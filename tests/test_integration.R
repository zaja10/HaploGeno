library(HaploGeno)
library(ggplot2)

message("--- Starting Integration Test ---")

# --- 1. Mock Genotype File ---
# Create a dummy FBM
n_ind <- 50
n_snp <- 100
geno_file <- tempfile(fileext = ".rds")
backing <- tempfile()

# Use bigstatsr to create dummy
if (requireNamespace("bigstatsr", quietly = TRUE)) {
    library(bigstatsr)
    X <- FBM(n_ind, n_snp, type = "double", backingfile = backing)
    X[] <- matrix(rnorm(n_ind * n_snp), n_ind, n_snp)
    sample_ids <- paste0("G", 1:n_ind)
    # We need to save the FBM correctly or ensure import handles it.
    # The import_genotypes expects .bed, .vcf, or .rds of bigSNP/FBM.
    # saveRDS(X) works for FBM.
    saveRDS(X, geno_file)
} else {
    stop("bigstatsr not installed.")
}

# Mock FA Model (structure mimic)
fa_mock <- list(
    fast = data.frame(
        Genotype = paste0("G", 1:n_ind),
        OP = rnorm(n_ind),
        RMSD = runif(n_ind)
    )
)
class(fa_mock) <- "fa_model"

# --- 2. Run Bridge ---
message("Testing init_from_met...")
# Case A: Full Match
haplo <- init_from_met(fa_mock, geno_path = geno_file, trait = "OP")
haplo$sample_ids <- paste0("G", 1:n_ind) # Mock sample IDs since FBM RDS doesn't store them efficiently without bigSNP wrapper usually

print(haplo)
# Expect: Phenotypes Loaded [x]
cat("Phenotype Length:", length(haplo$pheno), "\n")

if (length(haplo$pheno) != n_ind) stop("Phenotype loading failed.")

# --- 3. Verify Plotting (Mock Significance) ---
message("Testing Plotting...")
haplo$significance <- data.frame(
    BlockID = 1:10,
    P_Value = runif(10, 0, 0.1)
)
haplo$blocks <- data.frame(Start = 1:10, End = 1:10) # dummy

# S3 Wrapper
# Check if ggplot is returned
p <- plot(haplo, type = "manhattan")

if ("ggplot" %in% class(p)) {
    message("Plot returned ggplot object: SUCCESS")
} else {
    stop("Plot did not return ggplot object.")
}

message("--- Integration Test PASSED ---")
