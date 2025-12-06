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
  cat("========================================================\n")
  cat("           HAPLOGENO BREEDING REPORT                    \n")
  cat("========================================================\n\n")

  # 1. Overview
  n_samples <- if (!is.null(object$geno)) nrow(object$geno) else 0
  n_markers <- if (!is.null(object$geno)) ncol(object$geno) else 0
  n_blocks <- if (!is.null(object$blocks)) nrow(object$blocks) else 0

  cat("OVERVIEW:\n")
  cat("  Samples:    ", n_samples, "\n")
  cat("  Markers:    ", n_markers, "\n")
  cat("  Haplotypes: ", n_blocks, " defined blocks\n")
  if (!is.null(object$model_info)) {
    cat("  Model:      Ridge Regression (lambda =", round(object$model_info$lambda, 2), ")\n")
  }
  cat("\n")

  # 2. Driver Analysis
  cat("TOP DRIVER BLOCKS (Variance Explained):\n")
  if (!is.null(object$significance) && (!is.null(object$significance$PVE) || !is.null(object$significance$PVE_Adj))) {
    # Prefer adjusted PVE, then raw PVE, then Variance
    if (!is.null(object$significance$PVE_Adj)) {
      val <- object$significance$PVE_Adj
      lab <- "PVE(Adj)"
    } else if (!is.null(object$significance$PVE)) {
      val <- object$significance$PVE
      lab <- "PVE(Raw)"
    } else {
      val <- object$significance$Variance
      lab <- "Variance"
    }

    # Sort
    ord <- order(val, decreasing = TRUE)
    top_5 <- head(ord, 5)

    # Print table
    cat(sprintf("  %-10s %-15s %-12s %-10s\n", "BlockID", "Position", "Chrom", lab))
    cat("  --------------------------------------------------------\n")

    for (i in top_5) {
      bid <- object$significance$BlockID[i]
      # Get Pos/Chr from blocks/map
      start <- object$blocks$Start[i]
      end <- object$blocks$End[i]

      # Map info
      chrom <- "NA"
      pos_str <- paste0(start, "-", end) # fallback

      if (!is.null(object$map) && "chr" %in% names(object$map)) {
        # Use active markers map
        real_start <- if (!is.null(object$active_markers)) object$active_markers[start] else start
        chrom <- as.character(object$map$chr[real_start])

        if ("pos" %in% names(object$map)) {
          real_end <- if (!is.null(object$active_markers)) object$active_markers[end] else end
          p1 <- object$map$pos[real_start]
          p2 <- object$map$pos[real_end]
          pos_str <- paste0(round(p1 / 1e6, 2), "-", round(p2 / 1e6, 2), " Mb")
        }
      }

      v_print <- if (lab == "Variance") sprintf("%.4f", val[i]) else sprintf("%.2f%%", val[i] * 100)

      cat(sprintf("  %-10s %-15s %-12s %-10s\n", bid, pos_str, chrom, v_print))
    }

    # Total PVE
    total_pve <- sum(val, na.rm = TRUE)
    if (lab != "Variance") {
      cat("\n  Total Genomic Variance Explained: ", sprintf("%.2f%%", min(total_pve * 100, 100)), "\n")
    }
  } else {
    cat("  [Run calculate_pve() to see top drivers]\n")
  }
  cat("\n")

  # 3. Structural Analysis
  # We can't easily access private fa_results, so we check for public methods or just infer
  # If calculate_local_gebv is done, we have variance structure.

  cat("STATUS:\n")
  status_vec <- c(
    "Genotypes" = !is.null(object$geno),
    "Map" = !is.null(object$map),
    "Phenotypes" = !is.null(object$pheno),
    "Blocks" = !is.null(object$blocks),
    "Effects" = !is.null(object$marker_effects),
    "PVE" = !is.null(object$significance) && !is.null(object$significance$PVE)
  )

  # Print status in 2 columns
  for (i in seq(1, length(status_vec), by = 2)) {
    n1 <- names(status_vec)[i]
    s1 <- ifelse(status_vec[i], "[x]", "[ ]")

    if (i + 1 <= length(status_vec)) {
      n2 <- names(status_vec)[i + 1]
      s2 <- ifelse(status_vec[i + 1], "[x]", "[ ]")
      cat(sprintf("  %-20s %-20s\n", paste(s1, n1), paste(s2, n2)))
    } else {
      cat(sprintf("  %-20s\n", paste(s1, n1)))
    }
  }

  invisible(list(
    stats = list(samples = n_samples, markers = n_markers),
    status = status_vec
  ))
}
