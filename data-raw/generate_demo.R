devtools::load_all(".")
library(bigstatsr)

set.seed(42)

# Parameters
n <- 50
m <- 500
out_dir <- "inst/extdata"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Generating demo data...")

# Genotypes
X <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n, ncol = m)
# Add some structure/LD
for (i in seq(1, m - 10, by = 10)) {
    X[, i + 1] <- X[, i] # Perfect LD
}

# Map
map <- data.frame(
    chr = rep(1, m),
    id = paste0("snp", 1:m),
    pos = 1:m,
    ref = "A",
    alt = "G"
)

# Phenotypes (Signal in first 50 markers)
eff <- rnorm(50)
y <- as.vector(X[, 1:50] %*% eff + rnorm(n))

# Initialize
bk_file <- file.path(out_dir, "demo_data")
if (file.exists(paste0(bk_file, ".bk"))) file.remove(paste0(bk_file, ".bk"))
if (file.exists(paste0(bk_file, ".rds"))) file.remove(paste0(bk_file, ".rds"))

haplo <- HaploObject$new(bk_file)
haplo$import_genotypes(X)
haplo$load_map(map)
haplo$load_pheno(y)

# Process
message("Running pipeline...")
haplo$filter_monomorphic()
haplo$define_blocks_ld(r2_threshold = 0.5, window_size = 50)
haplo$encode_haplotypes()
haplo$compute_hrm()
haplo$estimate_marker_effects()
haplo$calculate_local_gebv()
haplo$test_significance()

# Save
message("Saving project...")
haplo$save_project(file.path(out_dir, "demo_data.rds"))

message("Done!")
