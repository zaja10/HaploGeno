#' Analyze Multi-Trait Haplotype Structure
#'
#' Performs Factor Analysis across multiple traits to identify pleiotropic genomic drivers.
#'
#' @param haplo_list A list of trained \code{HaploObject} instances (one per trait).
#' @param factors Number of latent factors to extract (default = 2).
#' @param top_n Number of top variance blocks to use per trait (default = 500).
#' @return A list containing Factor Analysis results (Loadings, Scores, Communality).
#' @importFrom stats sd varimax
#' @export
analyze_multitrait_structure <- function(haplo_list, factors = 2, top_n = 500) {
    # 1. Validation
    if (!is.list(haplo_list) || length(haplo_list) < 2) {
        stop("haplo_list must be a list of at least 2 HaploObjects.")
    }

    ref_obj <- haplo_list[[1]]
    if (!inherits(ref_obj, "HaploObject")) stop("Elements must be HaploObjects.")

    # Check consistency of samples
    ref_samples <- ref_obj$sample_ids

    for (i in 2:length(haplo_list)) {
        obj <- haplo_list[[i]]
        if (!inherits(obj, "HaploObject")) stop("Elements must be HaploObjects.")
        if (!identical(obj$sample_ids, ref_samples)) {
            stop("All HaploObjects must have identical Sample IDs.")
        }
    }

    message("Analyzing multi-trait structure across ", length(haplo_list), " traits...")

    # 2. Extract and Prepare Matrices
    # We want to select high-variance blocks from EACH trait.
    # It allows finding factors that might be driven by different blocks in different traits?
    # Or should we use the union of top blocks?
    # Let's simple cbind selected blocks.

    combined_matrix_list <- list()
    col_names <- character()
    block_indices <- integer()
    trait_indices <- integer()

    for (i in seq_along(haplo_list)) {
        obj <- haplo_list[[i]]
        if (is.null(obj$local_gebv)) stop("Local GEBVs not calculated for trait ", i)

        mat <- obj$local_gebv$matrix
        vars <- apply(mat, 2, var)

        # Select top N for this trait
        n_blocks <- length(vars)
        eff_n <- min(top_n, n_blocks)
        top_idx <- order(vars, decreasing = TRUE)[1:eff_n]

        sub_mat <- mat[, top_idx, drop = FALSE]

        # Store
        combined_matrix_list[[i]] <- sub_mat

        # Meta tracking
        trait_name <- paste0("Trait", i)
        col_names <- c(col_names, paste0(trait_name, "_B", top_idx))
        block_indices <- c(block_indices, top_idx)
        trait_indices <- c(trait_indices, rep(i, eff_n))
    }

    # Bind
    X_combined <- do.call(cbind, combined_matrix_list)
    colnames(X_combined) <- col_names

    # 3. Factor Analysis (SVD on Correlation/Scaled Matrix)
    # Scale columns
    X_scaled <- scale(X_combined)

    # Check for zero var / NA
    valid_cols <- which(attr(X_scaled, "scaled:scale") > 1e-10)
    if (length(valid_cols) < factors) stop("Not enough valid blocks for Factor Analysis.")

    X_scaled <- X_scaled[, valid_cols, drop = FALSE]

    # Filter meta info
    block_indices <- block_indices[valid_cols]
    trait_indices <- trait_indices[valid_cols]

    message("Factor Analysis on ", ncol(X_scaled), " combined block-traits...")

    # SVD
    s <- svd(X_scaled)

    U <- s$u[, 1:factors, drop = FALSE]
    D <- diag(s$d[1:factors], nrow = factors, ncol = factors)
    V <- s$v[, 1:factors, drop = FALSE]

    # Loadings
    scale_factor <- sqrt(nrow(X_scaled) - 1)
    L_unrotated <- V %*% D / scale_factor

    # Varimax
    vm <- varimax(L_unrotated)
    L_rotated <- vm$loadings
    class(L_rotated) <- "matrix" # strip 'loadings' class

    # Scores
    S_unrotated <- U * scale_factor
    S_rotated <- S_unrotated %*% vm$rotmat
    if (!is.null(ref_samples)) rownames(S_rotated) <- ref_samples

    colnames(L_rotated) <- paste0("Factor", 1:factors)
    colnames(S_rotated) <- paste0("Factor", 1:factors)

    # Communality
    h2 <- rowSums(L_rotated^2)
    u2 <- 1 - h2

    # Return structure
    structure(list(
        Loadings = L_rotated,
        Scores = S_rotated,
        Communality = h2,
        SpecificVar = u2,
        BlockIndices = block_indices,
        TraitIndices = trait_indices,
        Rotation = vm$rotmat
    ), class = "MultiTraitFA")
}

#' Print Multi-Trait FA Results
#' @param x A MultiTraitFA object.
#' @param ... Additional arguments.
#' @export
print.MultiTraitFA <- function(x, ...) {
    cat("Multi-Trait Factor Analysis Result\n")
    cat("Factors:", ncol(x$Scores), "\n")
    cat("Items (Block-Traits):", nrow(x$Loadings), "\n")
    cat("Traits Involved:", length(unique(x$TraitIndices)), "\n")
    invisible(x)
}
