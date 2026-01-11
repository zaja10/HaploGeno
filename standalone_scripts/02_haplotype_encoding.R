# Script: 02_haplotype_encoding.R
# Description: Encodes Study Genotypes into Haplotype IDs.
#              Includes an optimized C++ implementation and a pure R fallback.
# Dependencies: Rcpp, bigstatsr, future, future.apply

library(Rcpp)
library(future)
library(future.apply)
library(bigstatsr)

# --- 1. Inline C++ Function (Optimized) ---
# Compiles on source. Methodology: Order of First Appearance.
cpp_src <- "
#include <Rcpp.h>
#include <map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector encode_block_fast(NumericMatrix mat) {
    int n = mat.nrow();
    int m = mat.ncol();

    // Key: Haplotype (vector of ints), Value: ID (int)
    std::map<std::vector<int>, int> haplo_map;
    IntegerVector ids(n);

    int current_id = 1;

    for (int i = 0; i < n; i++) {
        // Construct key for current sample
        std::vector<int> row_key(m);
        for (int j = 0; j < m; j++) {
            row_key[j] = (int)mat(i, j);
        }

        // Check if haplotype exists
        std::map<std::vector<int>, int>::iterator it = haplo_map.lower_bound(row_key);
        
        if (it == haplo_map.end() || haplo_map.key_comp()(row_key, it->first)) {
            // New haplotype found: Assign next available ID
            haplo_map.insert(it, std::make_pair(row_key, current_id));
            ids[i] = current_id;
            current_id++;
        } else {
            // Existing haplotype: Retrieve ID
            ids[i] = it->second;
        }
    }
    return ids;
}
"

# Compile the C++ function immediately upon sourcing this script
message("Compiling C++ encoding function...")
sourceCpp(code = cpp_src)

# --- 2. Pure R Fallback Function ---
# Matches C++ logic: IDs assigned based on order of first appearance.
encode_block_R <- function(mat) {
    # 1. Collapse rows into string keys (e.g., "0_1_2")
    keys <- apply(mat, 1, paste, collapse = "_")
    
    # 2. Assign integer IDs
    # unique(keys) preserves the order of appearance.
    # match() finds the index in that unique vector.
    u_keys <- unique(keys)
    return(match(keys, u_keys))
}

# --- 3. Main Wrapper Function ---
#' Encode Haplotypes
#'
#' @param geno_mat Genotype matrix or FBM.
#' @param blocks_df Data frame of blocks (BlockID, Start, End).
#' @param n_cores Number of cores for parallel processing.
#' @param use_cpp Logical. Use C++ optimization (default: TRUE).
#' @return A matrix of Haplotype IDs (Rows = Samples, Cols = Blocks).
encode_haplotypes_func <- function(geno_mat, blocks_df, n_cores = 1, use_cpp = TRUE) {
    n_blocks <- nrow(blocks_df)
    
    # --- Parallel Logic Check ---
    # inline C++ functions (sourceCpp) do not export easily to parallel workers 
    # in 'multisession' (Windows) or 'future' environments without packaging.
    # To ensure stability, we force sequential mode if C++ is used.
    if (use_cpp && n_cores > 1) {
        message("Note: Parallel processing disabled for C++ mode to ensure stability.")
        message("      (C++ is sufficiently fast that sequential execution is recommended).")
        n_cores <- 1
    }

    message("Encoding haplotypes for ", n_blocks, " blocks (Cores: ", n_cores, ", Method: ", ifelse(use_cpp, "C++", "Pure R"), ")...")

    # Setup Plan
    if (n_cores > 1) {
        plan(multisession, workers = n_cores)
    } else {
        plan(sequential)
    }

    # Execute
    results <- future_lapply(1:n_blocks, function(i) {
        start_ind <- blocks_df$Start[i]
        end_ind <- blocks_df$End[i]

        # FBM vs Matrix handling is automatic with [ , ]
        # drop=FALSE ensures we pass a matrix even if block size is 1
        sub_mat <- geno_mat[, start_ind:end_ind, drop = FALSE]

        if (use_cpp) {
            return(encode_block_fast(sub_mat))
        } else {
            return(encode_block_R(sub_mat))
        }
    }, future.seed = TRUE)

    # Reset plan
    plan(sequential)

    # Combine columns
    haplo_mat <- do.call(cbind, results)
    colnames(haplo_mat) <- paste0("B", blocks_df$BlockID)
    
    message("Haplotype encoding complete.")
    return(haplo_mat)
}
