
#' Print method for HaploObject
#' @param x A HaploObject
#' @param ... Additional arguments
#' @export
print.HaploObject <- function(x, ...) {
  cat("HaploObject\n")
  cat("-----------\n")
  if (!is.null(x$geno)) {
    cat("Genotypes: ", nrow(x$geno), " samples x ", ncol(x$geno), " markers\n")
  } else {
    cat("Genotypes: Not loaded\n")
  }
  
  if (!is.null(x$blocks)) {
    cat("Blocks: ", nrow(x$blocks), " defined\n")
  } else {
    cat("Blocks: Not defined\n")
  }
  
  if (!is.null(x$hrm)) {
    cat("HRM: Computed\n")
  } else {
    cat("HRM: Not computed\n")
  }
  
  if (!is.null(x$significance)) {
    cat("Significance: Tested\n")
  }
}

#' Summary method for HaploObject
#' @param object A HaploObject
#' @param ... Additional arguments
#' @export
summary.HaploObject <- function(object, ...) {
  print(object)
  # Add more stats if needed
}
