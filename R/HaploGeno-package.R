#' @keywords internal
"_PACKAGE"

#' @useDynLib HaploGeno, .registration = TRUE
#' @import methods
#' @import ggplot2
#' @importFrom MASS ginv
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom Rcpp sourceCpp
NULL

# Silence R CMD check notes for Rcpp exported functions
utils::globalVariables(c("_HaploGeno_encode_block_fast", "_HaploGeno_ridge_solver_cpp"))
