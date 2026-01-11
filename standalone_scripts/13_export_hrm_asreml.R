# Script: 13_export_hrm_asreml.R
# Description: Converts Haplotype Matrix into ASReml-ready G-Inverse (.giv).
#              Calculates IBS, bends the matrix for inversion, and exports sparse format.
# Dependencies: Matrix, data.table

library(Matrix)
library(data.table)

#' Export Haplotype Relationship Matrix for ASReml
#'
#' @param haplo_mat Matrix of Haplotype IDs (Rows=Ind, Cols=Blocks). From Script 02.
#' @param out_prefix Filename prefix for output (e.g., "output/HaploG").
#' @param bending_factor Small value added to diagonal to ensure invertibility (default 0.001).
#' @return Invisible NULL. Saves .giv and .map files.
export_hrm_for_asreml <- function(haplo_mat, out_prefix = "HaploG", bending_factor = 0.001) {
    message("--- Generaring ASReml Haplotype Matrix ---")

    # 1. Compute IBS Relationship Matrix (Similarity)
    # Logic: Proportion of blocks where alleles are identical.
    # Note: For large matrices, this simple loop is safer than massive outer products.
    # If you have Script 07 loaded, you can use compute_hrm_func() instead.

    n_ind <- nrow(haplo_mat)
    n_blk <- ncol(haplo_mat)
    geno_names <- rownames(haplo_mat)

    if (is.null(geno_names)) stop("haplo_mat must have rownames (Genotype IDs).")

    message("Computing similarity for ", n_ind, " individuals...")

    # Efficient computation of IBS matrix (G)
    # G[i,j] = sum(haplo_mat[i,] == haplo_mat[j,]) / n_blk
    G <- matrix(0, nrow = n_ind, ncol = n_ind)
    rownames(G) <- colnames(G) <- geno_names

    # Using basic R for clarity, but could be C++ optimized
    # For user ease: efficient R approach
    for (i in seq_len(n_ind)) {
        for (j in i:n_ind) {
            # Count matches
            matches <- sum(haplo_mat[i, ] == haplo_mat[j, ], na.rm = TRUE)
            val <- matches / n_blk
            G[i, j] <- val
            G[j, i] <- val # Symmetric
        }
    }

    # 2. Bending (Ensure Invertibility)
    # HRMs are often singular because some lines share identical haplotypes.
    # We add a small noise to the diagonal (Ridge).
    message("Bending matrix (factor: ", bending_factor, ")...")
    G_bend <- G + diag(bending_factor, n_ind)

    # 3. Invert the Matrix
    message("Inverting G matrix...")
    # standard solve might fail if still singular, but bending handles most cases.
    G_inv <- solve(G_bend)

    # 4. Format for ASReml (.giv format)
    # ASReml needs a sparse matrix: RowIndex, ColIndex, Value
    # Only Lower Triangle

    # Create the Map (Genotype Names -> 1:N)
    map_df <- data.frame(
        ID = geno_names,
        Index = seq_len(n_ind),
        stringsAsFactors = FALSE
    )

    # Extract lower triangle of Inverse
    lt <- lower.tri(G_inv, diag = TRUE)
    g_sparse <- data.frame(
        Row = row(G_inv)[lt],
        Col = col(G_inv)[lt],
        Value = G_inv[lt]
    )

    # 5. Write Files
    giv_file <- paste0(out_prefix, ".giv")
    map_file <- paste0(out_prefix, ".map") # CSV format for reference

    message("Writing ", giv_file, " and ", map_file, "...")

    # Write .giv (Space separated, no headers usually preferred by ASReml, but check manual)
    write.table(g_sparse, file = giv_file, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)

    # Write Map
    write.csv(map_df, file = map_file, row.names = FALSE, quote = FALSE)

    message("Done. Use in ASReml with vm(Genotype, source='", giv_file, "')")
}
