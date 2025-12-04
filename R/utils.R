#' Load Demo Dataset
#'
#' Loads a pre-processed demo dataset for testing and visualization.
#' The dataset contains 50 samples and 500 markers, with blocks defined,
#' haplotypes encoded, and local GEBVs calculated.
#'
#' @return A HaploObject instance.
#' @export
load_demo_data <- function() {
    # Find files in package
    rds_file <- system.file("extdata", "demo_data.rds", package = "HaploGeno")
    bk_file <- system.file("extdata", "demo_data.bk", package = "HaploGeno")

    if (rds_file == "" || bk_file == "") {
        stop("Demo data not found. Please reinstall the package.")
    }

    # Copy to temp directory to allow write access and avoid messing with package install
    temp_dir <- tempdir()
    temp_rds <- file.path(temp_dir, "demo_data.rds")
    temp_bk <- file.path(temp_dir, "demo_data.bk")

    file.copy(rds_file, temp_rds, overwrite = TRUE)
    file.copy(bk_file, temp_bk, overwrite = TRUE)

    message("Loading demo data from temporary directory...")
    return(load_haplo_project(temp_rds))
}
