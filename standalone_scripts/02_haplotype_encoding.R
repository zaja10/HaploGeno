# Script: 02_haplotype_encoding.R
# Description: Encodes Study Genotypes into Haplotype IDs using C++ optimization.
# Dependencies: Rcpp, bigstatsr, future.apply

library(Rcpp)
library(future.apply)
library(bigstatsr)

# --- Inline C++ Function ---
# This removes the need for package compilation.
# We compile it on source.
cpp_src <- "
#include <Rcpp.h>
#include <map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector encode_block_fast(NumericMatrix mat) {
   int n = mat.nrow();
   int m = mat.ncol();

   // Use vector<int> as key to avoid string conversion overhead
   std::map<std::vector<int>, int> haplo_map;
   IntegerVector ids(n);

   int current_id = 1;

   for (int i = 0; i < n; i++) {
     // Construct vector key for current row
     std::vector<int> row_key(m);
     for (int j = 0; j < m; j++) {
       row_key[j] = (int)mat(i, j);
     }

     // Lookup or Insert
     std::map<std::vector<int>, int>::iterator it = haplo_map.lower_bound(row_key);
     if (it == haplo_map.end() || haplo_map.key_comp()(row_key, it->first)) {
       haplo_map.insert(it, std::make_pair(row_key, current_id));
       ids[i] = current_id;
       current_id++;
     } else {
       ids[i] = it->second;
     }
   }

   return ids;
}
"

# Compile C++ function
message("Compiling C++ encoding function...")
Rcpp::sourceCpp(code = cpp_src)

#' Encode Haplotypes
#'
#' @param geno_mat Genotype matrix or FBM.
#' @param blocks_df Data frame of blocks (BlockID, Start, End).
#' @param n_cores Number of cores for parallel processing (default: 1).
#' @return A matrix of Haplotype IDs (Rows = Samples, Cols = Blocks).
# --- Pure R Fallback Function ---
encode_block_R <- function(mat) {
    # Create a string key for each row (sample)
    # This is slower than C++ but uses only base R
    keys <- apply(mat, 1, paste, collapse = "_")
    # Convert to integer IDs based on unique keys
    # factor() assigns integers 1..K to levels
    # We want appearance order? factor uses alphabetical by default.
    # To match C++ (which assigns 1..K by appearance usually, or map order),
    # we should use match(keys, unique(keys)) to preserve appearance order?
    # The C++ map implementation:
    # std::map sorts keys. So IDs are assigned based on sorted key order?
    # Actually, the C++ code assigns 'current_id' (1, 2...) in order of APPEARANCE of unique haplotypes?
    # No, map inserts sorted.
    # Wait, let's check C++ logic:
    # it loop i=0..n.
    # if key not found, insert(key, current_id). current_id++.
    # So IDs are assigned in order of FIRST APPEARANCE in the dataset.
    # Correct pure R equivalent:
    u_keys <- unique(keys)
    return(match(keys, u_keys))
}

#' Encode Haplotypes
#'
#' @param geno_mat Genotype matrix or FBM.
#' @param blocks_df Data frame of blocks (BlockID, Start, End).
#' @param n_cores Number of cores for parallel processing (default: 1).
#' @param use_cpp Logical. Use C++ optimization (default: TRUE). If FALSE, uses pure R.
#' @return A matrix of Haplotype IDs (Rows = Samples, Cols = Blocks).
encode_haplotypes_func <- function(geno_mat, blocks_df, n_cores = 1, use_cpp = TRUE) {
    n_blocks <- nrow(blocks_df)
    message("Encoding haplotypes for ", n_blocks, " blocks (Parallel: ", n_cores > 1, ", Method: ", ifelse(use_cpp, "C++", "Pure R"), ")...")

    # Setup Parallel Plan
    if (n_cores > 1) {
        future::plan(future::multisession, workers = n_cores)
    } else {
        future::plan(future::sequential)
    }

    # Loop (Parallel)
    results <- future.apply::future_lapply(1:n_blocks, function(i) {
        start_ind <- blocks_df$Start[i]
        end_ind <- blocks_df$End[i]

        # Extract sub-matrix
        # handle FBM vs Matrix
        if (inherits(geno_mat, "FBM")) {
            sub_mat <- geno_mat[, start_ind:end_ind, drop = FALSE]
        } else {
            sub_mat <- geno_mat[, start_ind:end_ind, drop = FALSE]
        }

        if (use_cpp) {
            return(encode_block_fast(sub_mat))
        } else {
            return(encode_block_R(sub_mat))
        }
    }, future.seed = TRUE)

    # Combine results (Column bind)
    haplo_mat <- do.call(cbind, results)
    colnames(haplo_mat) <- sprintf("Block_%d", seq_len(ncol(haplo_mat)))

    message("Haplotype encoding complete.")
    return(haplo_mat)
}
