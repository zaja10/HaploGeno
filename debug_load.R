devtools::load_all(".")
library(bigstatsr)

temp_dir <- tempdir()
proj_name <- "debug_proj"
rds_path <- file.path(temp_dir, paste0(proj_name, ".rds"))
bk_path <- file.path(temp_dir, proj_name)

cat("Creating project...\n")
fbm <- FBM(10, 5, type = "double", backingfile = bk_path)
fbm[] <- 1

haplo <- HaploObject$new(bk_path)
haplo$geno <- fbm
haplo$save_project(rds_path)

cat("Moving project...\n")
new_dir <- file.path(temp_dir, "moved_debug")
if (!dir.exists(new_dir)) dir.create(new_dir)

file.rename(rds_path, file.path(new_dir, paste0(proj_name, ".rds")))
file.rename(paste0(bk_path, ".bk"), file.path(new_dir, paste0(proj_name, ".bk")))

new_rds_path <- file.path(new_dir, paste0(proj_name, ".rds"))
cat("Loading from: ", new_rds_path, "\n")

# Call internal function? No, it's exported.
loaded <- load_haplo_project(new_rds_path)

cat("Checking geno...\n")
print(loaded$geno$backingfile)
tryCatch(
    {
        print(loaded$geno[1, 1])
        cat("Success!\n")
    },
    error = function(e) {
        cat("Error: ", e$message, "\n")
    }
)
