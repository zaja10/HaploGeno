#' Initialize HaploObject from UtilityFunctions Output
#'
#' Creates a HaploObject using metrics (OP, RMSD, Index) from a phenotypic analysis
#' as the target phenotype. Performs automatic ID alignment.
#'
#' @param met_object An object of class \code{fa_model} (single trait) or
#'                   \code{mt_selection_results} (multi-trait index).
#' @param geno_path Path to the genotype file (VCF/BED/RDS).
#' @param trait Character. If input is \code{fa_model}, specifying "OP", "RMSD", or "Factor".
#'        If input is \code{mt_selection_results}, specify "Index" or a specific trait name.
#' @param metric_col Character (Optional). Manually specify column name if auto-detection fails.
#' @param factor_idx Integer. If trait="Factor", which factor index to use? Default 1.
#'
#' @return A \code{HaploObject} with genotypes mapped and phenotypes loaded.
#' @export
init_from_met <- function(met_object, geno_path, trait = "Index", metric_col = NULL, factor_idx = 1) {
    # 1. Extract Phenotype Dataframe
    df <- NULL
    target_col <- NULL

    if (inherits(met_object, "fa_model")) {
        if (trait %in% c("OP", "RMSD")) {
            if (is.null(met_object$fast)) stop("Run calculate_fast_indices() on the FA model first.")
            df <- met_object$fast
            target_col <- trait
        } else if (trait == "Factor") {
            # Extract factor scores if needed
            if (is.null(met_object$scores$rotated)) stop("No rotated scores found in FA model.")
            df <- as.data.frame(met_object$scores$rotated)
            # Scores matrix usually has Genotypes as rownames
            df$Genotype <- rownames(df)

            # Name of factor column in the matrix?
            # Actually scores$rotated is a matrix probably.
            # Let's check dimensions
            k <- ncol(met_object$scores$rotated)
            if (factor_idx > k) stop("Requested factor index exceeds number of factors.")

            # The matrix columns might be named "Fac1", "Fac2" etc or just indices
            target_col <- colnames(met_object$scores$rotated)[factor_idx]
            if (is.null(target_col)) target_col <- paste0("Fac", factor_idx)

            # Ensure the dataframe has this column
            # If it was a matrix, as.data.frame should keep colnames
            # If not, rename
            if (!target_col %in% names(df)) {
                names(df)[factor_idx] <- target_col
            }
        } else {
            # Maybe the user passed a trait name but meant OP/RMSD?
            # Or passed "Index" but provided an fa_model?
            if (trait == "Index") stop("fa_model input does not contain 'Index'. Run rank_genotypes() or use calculate_mt_index().")
            stop("For fa_model, trait must be 'OP', 'RMSD', or 'Factor'.")
        }
    } else if (inherits(met_object, "mt_selection_results") || inherits(met_object, "data.frame")) {
        # Support raw dataframe if it has Genotype column
        if (inherits(met_object, "mt_selection_results")) df <- met_object$index else df <- met_object

        if (!"Genotype" %in% names(df)) stop("Input dataframe must contain a 'Genotype' column.")

        if (trait == "Index") {
            target_col <- "Index"
        } else {
            # Try to find specific trait column (e.g. Yield_OP_Z)
            # If user passed "Yield", look for it or Yield_OP_Z?
            if (trait %in% names(df)) {
                target_col <- trait
            } else {
                # Try constructing typical names
                candidates <- c(paste0(trait, "_OP_Z"), paste0(trait, "_Comp"))
                found <- intersect(candidates, names(df))
                if (length(found) > 0) target_col <- found[1]
            }
        }
    } else {
        stop("Unknown object class. Must be 'fa_model', 'mt_selection_results', or 'data.frame'.")
    }

    if (!is.null(metric_col)) target_col <- metric_col

    if (is.null(target_col) || is.null(df[[target_col]])) {
        stop(paste("Could not find trait/column:", trait, "(Target:", target_col, ")"))
    }

    message(paste(">>> Importing Phenotypes:", target_col))

    # 2. Initialize HaploObject
    haplo <- HaploObject$new(tempfile())
    haplo$import_genotypes(geno_path)

    # 3. ID Alignment (The Critical Step)
    geno_ids <- haplo$sample_ids
    pheno_ids <- df$Genotype

    # Ensure unique phenotypes (if multiple rows per genotype, stop or average?)
    if (any(duplicated(pheno_ids))) {
        warning("Duplicate genotypes found in phenotype data. Taking the first occurrence.")
        df <- df[!duplicated(df$Genotype), ]
        pheno_ids <- df$Genotype
    }

    pheno_vals <- setNames(df[[target_col]], pheno_ids)

    # Intersection
    common <- intersect(geno_ids, pheno_ids)
    if (length(common) == 0) stop("No matching IDs found between Genotypes and Phenotypes!")

    coverage <- (length(common) / length(geno_ids)) * 100
    message(sprintf("Matched %d genotypes (%.1f%% of genomic file).", length(common), coverage))

    if (coverage < 10) warning("Very low ID alignment (<10%). Check ID formats (e.g. spaces, capitalization).")

    # Align Phenotype Vector to Genotype Rows (Inserting NAs for missing lines)
    # This ensures row 1 of geno corresponds to element 1 of pheno
    aligned_y <- pheno_vals[geno_ids]

    # Check for NAs in target set (genotypes present in file but missing pheno)
    n_missing <- sum(is.na(aligned_y))
    if (n_missing > 0) {
        message(sprintf("Note: %d genotypes in file have no phenotype (set to NA).", n_missing))
    }

    # 4. Filter missing phenotypes??
    # HaploGeno handles NAs?
    # Usually we want to SUBSET the genotypes to only those with phenotypes to save memory/time.
    # But FBM is immutable-ish. We can use 'ind.row' in bigstatsr functions.
    # HaploObject currently expects pheno to match geno.
    # If we want to support subsetting, we might need a subset_genotypes method or just load all.
    # For now, we load all.

    # 5. Load
    haplo$load_pheno(as.vector(aligned_y))

    return(haplo)
}
