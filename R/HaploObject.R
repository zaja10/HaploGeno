#' HaploObject Class
#'
#' @import R6
#' @import bigstatsr
#' @import data.table
#' @import future
#' @import future.apply
#' @import bigsnpr
#' @import progressr
#' @export
HaploObject <- R6::R6Class("HaploObject",
    public = list(
        #' @field geno Filebacked Big Matrix of genotypes.
        geno = NULL,
        #' @field map Data.table containing marker map information.
        map = NULL,
        #' @field pheno Vector of phenotypes.
        pheno = NULL,
        #' @field sample_ids Vector of sample IDs.
        sample_ids = NULL,
        #' @field blocks Data.table defining haplotype blocks.
        blocks = NULL,
        #' @field haplo_mat Matrix of encoded haplotypes.
        haplo_mat = NULL,
        #' @field hrm Haplotype Relationship Matrix.
        hrm = NULL,
        #' @field marker_effects Vector of estimated marker effects.
        marker_effects = NULL,
        #' @field local_gebv List containing local GEBV matrix and variances.
        local_gebv = NULL,
        #' @field significance Data.table of significance test results.
        significance = NULL,
        #' @field active_markers Integer vector of indices of markers with non-zero variance.
        active_markers = NULL,
        #' @description
        #' Initialize a new HaploObject.
        #' @param backing_file_path Path to the backing file for the FBM.
        initialize = function(backing_file_path) {
            # Initialize with a path for the FBM
            private$backing_file <- backing_file_path
        },

        #' @description
        #' Print a summary of the HaploObject.
        print = function() {
            cat("=== HaploObject Summary ===\n")

            # Samples
            if (!is.null(self$geno)) {
                cat(sprintf("Samples: %d\n", nrow(self$geno)))
                cat(sprintf("Markers: %d\n", ncol(self$geno)))
                if (!is.null(self$active_markers)) {
                    cat(sprintf(
                        "Active Markers: %d (%.1f%%)\n",
                        length(self$active_markers),
                        100 * length(self$active_markers) / ncol(self$geno)
                    ))
                }
            } else {
                cat("Samples: 0 (Genotypes not loaded)\n")
            }

            # Status
            cat("Status:\n")
            cat(sprintf("  [%s] Genotypes Loaded\n", ifelse(!is.null(self$geno), "x", " ")))
            cat(sprintf("  [%s] Map Loaded\n", ifelse(!is.null(self$map), "x", " ")))
            cat(sprintf("  [%s] Phenotypes Loaded\n", ifelse(!is.null(self$pheno), "x", " ")))
            cat(sprintf("  [%s] Blocks Defined\n", ifelse(!is.null(self$blocks), "x", " ")))
            cat(sprintf("  [%s] Haplotypes Encoded\n", ifelse(!is.null(self$haplo_mat), "x", " ")))
            cat(sprintf("  [%s] Model Trained\n", ifelse(!is.null(self$marker_effects), "x", " ")))

            invisible(self)
        },
        #' @description
        #' Import genotypes from a file or matrix.
        #' @param matrix_or_path A matrix or path to a file (bed/vcf/rds).
        import_genotypes = function(matrix_or_path) {
            # Logic to convert in-memory matrix to FBM or load from file
            if (is.matrix(matrix_or_path)) {
                # Check if backing file exists, if so, append or overwrite?
                # Currently, we default to creating a new backing file for each session to ensure data integrity.
                # Convert to FBM
                # We enforce type = "double" to ensure compatibility with downstream linear algebra operations (big_cor, big_scale).
                if (!is.null(rownames(matrix_or_path))) {
                    self$sample_ids <- rownames(matrix_or_path)
                }
                self$geno <- bigstatsr::as_FBM(matrix_or_path,
                    type = "double",
                    backingfile = private$backing_file
                )
            } else if (is.character(matrix_or_path)) {
                if (!file.exists(matrix_or_path)) stop("File not found: ", matrix_or_path)

                ext <- tools::file_ext(matrix_or_path)

                if (ext == "bed") {
                    # PLINK bed file
                    message("Importing PLINK .bed file...")
                    # snp_readBed creates a .rds file
                    rds_file <- bigsnpr::snp_readBed(matrix_or_path, backingfile = private$backing_file)
                    obj <- bigsnpr::snp_attach(rds_file)
                    self$geno <- obj$genotypes
                    # Also load map/fam if possible?
                    # bigsnpr object contains map/fam information.
                    # TODO: Extract map and store it automatically for better UX.
                    # For now, just getting genotypes.
                } else if (ext == "vcf" || ext == "gz") {
                    # VCF file (assuming .vcf or .vcf.gz)
                    message("Importing VCF file (this may take a while)...")
                    # snp_readVCF creates a .rds file
                    # Note: snp_readVCF might require specific formatting
                    rds_file <- bigsnpr::snp_readVCF(matrix_or_path, backingfile = private$backing_file)
                    obj <- bigsnpr::snp_attach(rds_file)
                    self$geno <- obj$genotypes
                } else if (ext == "csv" || ext == "txt") {
                    # CSV/TXT file with potential variable encodings
                    private$read_csv_genotypes(matrix_or_path)
                } else if (ext == "rds") {
                    # Bigstatsr/Bigsnpr RDS
                    obj <- bigsnpr::snp_attach(matrix_or_path)
                    if (inherits(obj, "bigSNP")) {
                        self$geno <- obj$genotypes
                        if (!is.null(obj$fam$sample.ID)) self$sample_ids <- obj$fam$sample.ID
                    } else {
                        if (file.exists(paste0(matrix_or_path, ".bk"))) {
                            self$geno <- bigstatsr::big_attach(paste0(matrix_or_path, ".rds"))
                        } else {
                            stop("Unsupported file format or backing file not found.")
                        }
                    }
                } else {
                    stop("File not found: ", matrix_or_path)
                }
            } else {
                stop("Invalid input for genotypes. Must be a matrix or a file path.")
            }
        },
        #' @description
        #' Load marker map.
        #' @param map_data A data.frame or path to a file.
        load_map = function(map_data) {
            # map_data can be a data.frame or a path to a file
            if (is.character(map_data) && file.exists(map_data)) {
                self$map <- data.table::fread(map_data)
            } else if (is.data.frame(map_data)) {
                self$map <- data.table::as.data.table(map_data)
            } else {
                stop("Invalid map data. Must be a data.frame or file path.")
            }

            # Validation: check number of markers matches genotypes
            if (!is.null(self$geno)) {
                if (nrow(self$map) != ncol(self$geno)) {
                    warning(
                        "Number of markers in map (", nrow(self$map),
                        ") does not match number of columns in genotypes (", ncol(self$geno), ")."
                    )
                }
            }

            # Check required columns
            required_cols <- c("chr", "id", "pos", "ref", "alt")
            missing_cols <- setdiff(required_cols, names(self$map))
            if (length(missing_cols) > 0) {
                warning(
                    "Map is missing standard columns: ", paste(missing_cols, collapse = ", "),
                    ". Some functions (e.g. Manhattan plot) may not work correctly."
                )
            }

            if ("pos" %in% names(self$map) && !is.numeric(self$map$pos)) {
                warning("'pos' column in map is not numeric. Attempting conversion.")
                self$map$pos <- as.numeric(self$map$pos)
            }
        },
        #' @description
        #' Load phenotypes.
        #' @param pheno_data A vector, data.frame, or path to a file.
        load_pheno = function(pheno_data) {
            # pheno_data can be a vector, data.frame, or path
            if (is.character(pheno_data) && file.exists(pheno_data)) {
                # Assume single column or specific format; simple fread for now
                dt <- data.table::fread(pheno_data)
                if (ncol(dt) == 1) {
                    self$pheno <- dt[[1]]
                } else {
                    # If multiple columns, assume 'y' or last column is phenotype
                    # This is a simplification; might need explicit col name arg
                    self$pheno <- dt[[ncol(dt)]]
                }
            } else if (is.vector(pheno_data)) {
                self$pheno <- pheno_data
            } else if (is.data.frame(pheno_data)) {
                # Assume last column
                self$pheno <- pheno_data[[ncol(pheno_data)]]
            } else {
                stop("Invalid phenotype data.")
            }

            # Validation
            if (!is.null(self$geno)) {
                if (length(self$pheno) != nrow(self$geno)) {
                    stop(
                        "Number of phenotypes (", length(self$pheno),
                        ") does not match number of samples in genotypes (", nrow(self$geno), ")."
                    )
                }
            }
        },
        #' @description
        #' Get a subset of the genotype matrix.
        #' @param row_ind Row indices.
        #' @param col_ind Column indices.
        get_subset = function(row_ind, col_ind) {
            return(self$geno[row_ind, col_ind])
        },
        #' @description
        #' Define haplotype blocks using fixed window size.
        #' @param window_size Number of markers per block.
        define_blocks_fixed = function(window_size) {
            # window_size in number of markers (e.g., 100)
            # Populates self$blocks with Start/End indices
            if (is.null(self$map)) stop("Map not loaded. Please load map first.")


            if (!is.null(self$active_markers)) {
                n_markers <- length(self$active_markers)
            } else {
                n_markers <- ncol(self$geno)
            }

            starts <- seq(1, n_markers, by = window_size)
            ends <- pmin(starts + window_size - 1, n_markers)

            self$blocks <- data.table::data.table(
                BlockID = seq_along(starts),
                Start = starts,
                End = ends
            )

            message("Defined ", nrow(self$blocks), " blocks.")
        },
        #' @description
        #' Define haplotype blocks using LD scan.
        #' @param r2_threshold r2 threshold for block definition.
        #' @param tolerance Number of failures allowed before ending a block.
        #' @param window_size Maximum window size to scan.
        define_blocks_ld = function(r2_threshold = 0.1, tolerance = 3, window_size = 2000) {
            if (is.null(self$geno)) stop("Genotypes not loaded.")

            message("Defining blocks using LD scan (r2 >= ", r2_threshold, ", tol = ", tolerance, ") with window size ", window_size, "...")

            message("Defining blocks using LD scan (r2 >= ", r2_threshold, ", tol = ", tolerance, ") with window size ", window_size, "...")

            if (!is.null(self$active_markers)) {
                n_markers <- length(self$active_markers)
            } else {
                n_markers <- ncol(self$geno)
            }
            starts <- c()
            ends <- c()

            current_start <- 1

            progressr::with_progress({
                p <- progressr::progressor(steps = n_markers)

                while (current_start <= n_markers) {
                    # Define window indices
                    window_end <- min(current_start + window_size, n_markers)
                    indices <- current_start:window_end

                    if (length(indices) < 2) {
                        # End of genome
                        starts <- c(starts, current_start)
                        ends <- c(ends, window_end)
                        p(amount = length(indices))
                        break
                    }

                    # Use big_crossprodSelf to compute X'X on scaled data
                    # This gives the correlation matrix * (n-1)
                    # Use big_crossprodSelf to compute X'X on scaled data
                    # This gives the correlation matrix * (n-1)

                    real_indices <- indices
                    if (!is.null(self$active_markers)) {
                        real_indices <- self$active_markers[indices]
                    }

                    K <- bigstatsr::big_crossprodSelf(self$geno, fun.scaling = bigstatsr::big_scale(), ind.col = real_indices)

                    # Convert to standard matrix and normalize
                    corr_mat <- K[] / (nrow(self$geno) - 1)

                    # Extract correlations with the first marker in the window (current_start)
                    r_vals <- corr_mat[1, ]
                    r2_vals <- r_vals^2

                    # Scan for block end
                    failures <- 0
                    block_end_rel <- 1 # Relative index in r2_vals

                    for (j in 2:length(r2_vals)) {
                        if (is.na(r2_vals[j])) r2_vals[j] <- 0

                        if (r2_vals[j] >= r2_threshold) {
                            block_end_rel <- j
                            failures <- 0
                        } else {
                            failures <- failures + 1
                        }

                        if (failures > tolerance) {
                            break
                        }
                    }

                    # Absolute end index
                    current_end <- indices[block_end_rel]

                    starts <- c(starts, current_start)
                    ends <- c(ends, current_end)

                    # Update progress
                    p(amount = current_end - current_start + 1)

                    current_start <- current_end + 1
                }
            })

            self$blocks <- data.table::data.table(
                BlockID = seq_along(starts),
                Start = starts,
                End = ends
            )

            message("Defined ", nrow(self$blocks), " blocks.")
        },
        #' @description
        #' Estimate marker effects using Ridge Regression.
        #' @param lambda Regularization parameter.
        estimate_marker_effects = function(lambda = 1.0) {
            if (is.null(self$geno)) stop("Genotypes not loaded.")
            if (is.null(self$pheno)) stop("Phenotypes not loaded.")

            message("Estimating marker effects using Ridge Regression...")

            # Prepare data
            # For large datasets, we should use big_spLinReg or similar, but that's Lasso/ElasticNet.
            # For Ridge, we can use big_univLinReg (GWAS) but that's univariate.
            # We need multivariate Ridge.

            # If N < P, use dual form: u = Z' (ZZ' + lambda I)^-1 y
            # If N > P, use primal form: u = (Z'Z + lambda I)^-1 Z'y

            n <- nrow(self$geno)
            p <- ncol(self$geno)
            y <- self$pheno

            # Center genotypes? Usually yes for genomic prediction.
            # We can use bigstatsr::big_scale() implicitly if we use bigstatsr functions.
            # But for manual algebra, we need to be careful.

            # Let's assume N < P (typical for genomic prediction)
            # We need K = ZZ' (Genomic Relationship Matrix)

            # Using bigstatsr to compute GRM efficiently
            # K = big_crossprodSelf(geno, fun.scaling = big_scale())
            # This returns K = X'X if X is n x p? No, crossprodSelf is X'X.
            # tcrossprodSelf is XX'.

            # big_tcrossprodSelf computes XX'
            # big_tcrossprodSelf computes XX'
            ind_col <- if (!is.null(self$active_markers)) self$active_markers else bigstatsr::cols_along(self$geno)

            # Robustness: Check for zero variance markers in the selection
            # big_scale() will error if any column has 0 variance.
            # Efficiently check variances using big_colstats
            stats_check <- bigstatsr::big_colstats(self$geno, ind.col = ind_col)
            keep_bool <- stats_check$var > 1e-8

            if (!all(keep_bool)) {
                n_drop <- sum(!keep_bool)
                message("Automatically excluding ", n_drop, " monomorphic markers from analysis.")
                ind_col <- ind_col[keep_bool]

                # Update active markers to persist this cleanup
                self$active_markers <- ind_col

                # If no markers left?
                if (length(ind_col) == 0) stop("No polymorphic markers left for analysis.")
            }

            K <- bigstatsr::big_tcrossprodSelf(self$geno, fun.scaling = bigstatsr::big_scale(), ind.col = ind_col)
            K <- K[] # Convert to standard matrix

            # Solve dual: alpha = (K + lambda I)^-1 y
            I <- diag(n)
            alpha <- solve(K + lambda * I, y)

            # Back-calculate marker effects: u = Z' alpha
            # u = Z_scaled' * alpha
            # Z_scaled = (Z - center) / scale
            # u = (Z' alpha - center * sum(alpha)) / scale

            stats <- bigstatsr::big_scale()(self$geno, ind.col = ind_col)
            center <- stats$center
            scale <- stats$scale

            raw_u <- bigstatsr::big_cprodVec(self$geno, alpha, ind.col = ind_col)
            sum_alpha <- sum(alpha)

            u_hat <- (raw_u - center * sum_alpha) / scale

            self$marker_effects <- as.vector(u_hat)
            message("Marker effects estimated.")
        },
        #' @description
        #' Encode haplotypes into integer IDs.
        #' @param n_cores Number of cores to use.
        encode_haplotypes = function(n_cores = 1) {
            if (is.null(self$blocks)) stop("Blocks not defined. Run define_blocks_* first.")
            if (is.null(self$geno)) stop("Genotypes not loaded.")

            n_blocks <- nrow(self$blocks)
            message("Encoding haplotypes for ", n_blocks, " blocks (Parallel: ", n_cores > 1, ")...")

            if (n_cores > 1) {
                future::plan(future::multisession, workers = n_cores)
            } else {
                future::plan(future::sequential)
            }

            progressr::with_progress({
                p <- progressr::progressor(steps = n_blocks)

                results <- future.apply::future_lapply(1:n_blocks, function(i) {
                    p()
                    start_ind <- self$blocks$Start[i]
                    end_ind <- self$blocks$End[i]

                    # FIX: Add drop = FALSE to prevent R from converting 1-column matrices to vectors
                    real_indices <- if (!is.null(self$active_markers)) self$active_markers[start_ind:end_ind] else start_ind:end_ind
                    sub_mat <- self$geno[, real_indices, drop = FALSE]

                    # Call the optimized C++ function
                    return(HaploGeno::encode_block_fast(sub_mat))
                }, future.seed = TRUE)
            })

            self$haplo_mat <- do.call(cbind, results)
            message("Haplotype encoding complete.")
        },
        #' @description
        #' Calculate local GEBVs.
        #' @param n_cores Number of cores to use.
        calculate_local_gebv = function(n_cores = 1) {
            if (is.null(self$blocks)) stop("Blocks not defined.")
            if (is.null(self$marker_effects)) stop("Marker effects not estimated.")

            message("Calculating localGEBV and variances (Parallel: ", n_cores > 1, ")...")

            n_blocks <- nrow(self$blocks)

            # --- OPTIMIZATION: Pre-calculate scaling factors globally ---
            message("Pre-calculating global scaling statistics...")
            # This scans the whole matrix ONCE
            ind_col <- if (!is.null(self$active_markers)) self$active_markers else bigstatsr::cols_along(self$geno)
            stats_all <- bigstatsr::big_scale()(self$geno, ind.col = ind_col)
            centers_all <- stats_all$center
            scales_all <- stats_all$scale

            if (n_cores > 1) {
                future::plan(future::multisession, workers = n_cores)
            } else {
                future::plan(future::sequential)
            }

            progressr::with_progress({
                p <- progressr::progressor(steps = n_blocks)

                results <- future.apply::future_lapply(1:n_blocks, function(i) {
                    p()
                    start <- self$blocks$Start[i]
                    end <- self$blocks$End[i]
                    indices <- start:end

                    # Handle active_markers (subsetting)
                    if (!is.null(self$active_markers)) {
                        # Logic Update: Treat start/end as relative indices into active_markers
                        # indices (start:end) are relative to the active_markers vector

                        # Indices in active_markers vector
                        rel_indices <- indices

                        # Map to absolute columns in genotype matrix
                        abs_indices <- self$active_markers[rel_indices]

                        # Marker effects and stats are indexed by relative position (1..N_active)
                        u_block <- self$marker_effects[rel_indices]
                        mus <- centers_all[rel_indices]
                        sigmas <- scales_all[rel_indices]

                        # big_prodVec needs absolute column indices
                        real_indices <- abs_indices
                    } else {
                        # All markers active
                        u_block <- self$marker_effects[indices]
                        mus <- centers_all[indices]
                        sigmas <- scales_all[indices]
                        real_indices <- indices
                    }

                    weights <- u_block / sigmas
                    offset <- sum(mus * weights)

                    # big_prodVec is efficient as it uses memory mapping
                    gebv <- bigstatsr::big_prodVec(self$geno, weights, ind.col = real_indices)
                    gebv <- gebv - offset

                    return(list(gebv = gebv, var = var(gebv)))
                }, future.seed = TRUE)
            })

            local_gebvs <- do.call(cbind, lapply(results, function(x) x$gebv))
            variances <- unlist(lapply(results, function(x) x$var))

            if (!is.null(self$sample_ids)) {
                rownames(local_gebvs) <- self$sample_ids
            }

            self$local_gebv <- list(
                matrix = local_gebvs,
                variances = variances
            )

            message("Local GEBV calculation complete.")
        },
        #' @description
        #' Get marker names for a specific block.
        #' @param block_id Integer ID of the block.
        get_block_markers = function(block_id) {
            if (is.null(self$blocks)) stop("Blocks not defined.")
            if (block_id < 1 || block_id > nrow(self$blocks)) stop("Invalid block ID.")

            start <- self$blocks$Start[block_id]
            end <- self$blocks$End[block_id]
            indices <- start:end

            if (!is.null(self$active_markers)) {
                abs_indices <- self$active_markers[indices]
            } else {
                abs_indices <- indices
            }

            if (!is.null(self$map) && "id" %in% names(self$map)) {
                return(self$map$id[abs_indices])
            } else {
                return(paste0("Marker_", abs_indices))
            }
        },
        #' @description
        #' Test significance of local GEBVs.
        #' @return Data.table of p-values.
        test_significance = function() {
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated.")

            message("Testing significance (Scaled Inverse Chi-Squared)...")

            vars <- self$local_gebv$variances

            # Scaled Inverse Chi-Squared Test
            # H0: variance = 0? Or variance = expected?
            # Shaffer et al. usually compares to a null distribution or uses a specific parameterization.
            # The prompt says: "Uses Scaled Inverse Chi-Squared (v=0.5, tau^2=0.5)"

            # This implies we are calculating a p-value based on the observed variance
            # under a Scaled Inv-Chi-Sq distribution with parameters nu=0.5, tau2=0.5.

            # The PDF of Scaled Inv-Chi-Sq(nu, tau2) is ...
            # Actually, usually we use this for Bayesian priors.
            # If used for testing, maybe we are testing if the variance is consistent with noise?
            # Or is it a p-value FROM the CDF?

            # Let's assume we want P(X > var) where X ~ Scaled-Inv-Chi-Sq(0.5, 0.5)
            # Scaled-Inv-Chi-Sq(nu, tau2) is equivalent to:
            # (nu * tau2) / X ~ Chi-Sq(nu)

            # So X ~ (nu * tau2) / Chi-Sq(nu)
            # We want P(X > observed_var)
            # = P( (nu * tau2) / Chi-Sq(nu) > observed_var )
            # = P( Chi-Sq(nu) < (nu * tau2) / observed_var )

            nu <- 0.5
            tau2 <- 0.5
            scale_param <- nu * tau2

            # Calculate statistic for each block
            # Avoid division by zero
            safe_vars <- vars
            safe_vars[safe_vars < 1e-10] <- 1e-10

            chi_sq_stats <- scale_param / safe_vars

            # P-value = P(Chi-Sq(nu) < stat)
            # This seems inverted. Usually large variance is significant.
            # If variance is LARGE, then stat is SMALL.
            # P(Chi-Sq < small) is small.
            # So yes, this direction makes sense: Small p-value for large variance.

            p_values <- pchisq(chi_sq_stats, df = nu)

            self$significance <- data.table::data.table(
                BlockID = self$blocks$BlockID,
                Variance = vars,
                P_Value = p_values
            )

            return(self$significance)
        },
        #' @description
        #' Calculate Percent Variance Explained (PVE) by each block.
        #' Computes r^2 between local GEBVs and phenotypes.
        #' @return Updated significance table.
        calculate_pve = function() {
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated.")
            if (is.null(self$pheno)) stop("Phenotypes not loaded.")

            message("Calculating Phenotypic Variance Explained (PVE)...")

            # cor() is fast and base R, handle missing values with use="pairwise"
            # local_gebv$matrix is N x M, pheno is N x 1
            # We want correlation of each block (column) with pheno

            # Ensure pheno is vector
            y <- as.vector(self$pheno)

            # Only compute for available data
            valid_idx <- which(!is.na(y))
            if (length(valid_idx) < length(y)) {
                warning("Phenotypes contain NAs. Using pairwise deletion for correlation.")
                y <- y[valid_idx]
                X <- self$local_gebv$matrix[valid_idx, , drop = FALSE]
            } else {
                X <- self$local_gebv$matrix
            }

            cor_vals <- cor(X, y)
            r2 <- as.vector(cor_vals^2)

            # Update significance table
            if (is.null(self$significance)) {
                # If significance table doesn't exist, create partial one
                self$test_significance()
            }

            self$significance$PVE <- r2
            return(self$significance)
        },
        #' @description
        #' Analyze block structure using internal Factor Analysis (Latent Regression).
        #' Models Local GEBVs as function of latent genomic gradients.
        #' @param top_n Number of top variance blocks to use.
        #' @param factors Number of latent factors to extract.
        analyze_block_structure = function(top_n = 500, factors = 2) {
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated.")

            message("Running Haplo-FA (Internal Factor Analysis) on top ", top_n, " blocks...")

            # 1. Select top blocks by variance
            # If significance table exists, use that. Else calculate variances.
            if (!is.null(self$significance)) {
                vars <- self$significance$Variance
            } else {
                vars <- apply(self$local_gebv$matrix, 2, var)
            }

            # Handle NAs in variance (shouldn't happen but be safe)
            vars[is.na(vars)] <- 0

            n_blocks <- length(vars)
            effective_n <- min(top_n, n_blocks)

            # Indices of top blocks
            top_indices <- order(vars, decreasing = TRUE)[1:effective_n]

            # Subset GEBV matrix
            X <- self$local_gebv$matrix[, top_indices, drop = FALSE]

            # 2. Scale (Z-score)
            X_scaled <- scale(X)
            # Remove columns with zero variance (if any slipped through)
            valid_cols <- which(attr(X_scaled, "scaled:scale") > 1e-10)
            if (length(valid_cols) < factors) stop("Not enough valid blocks for Factor Analysis.")

            X_scaled <- X_scaled[, valid_cols, drop = FALSE]
            final_indices <- top_indices[valid_cols]

            # 3. SVD for Factor Extraction
            # We want Factor Analysis, can approximate with PCA (SVD on correlation matrix)
            # Or SVD on data matrix directly.
            # X_scaled = U D V'.
            # Loadings = V * D / sqrt(N-1). Scores = U * sqrt(N-1)

            s <- svd(X_scaled)

            # Limit to 'factors' dimensions
            U <- s$u[, 1:factors, drop = FALSE]
            D <- diag(s$d[1:factors], nrow = factors, ncol = factors)
            V <- s$v[, 1:factors, drop = FALSE]

            # Initial (Unrotated) Loadings and Scores
            # Loadings: Correlation between Item (Block) and Factor
            # L = V %*% D / sqrt(nrow(X)-1)
            scale_factor <- sqrt(nrow(X_scaled) - 1)
            L_unrotated <- V %*% D / scale_factor

            # 4. Varimax Rotation
            # Improves interpretability (simple structure)
            vm <- varimax(L_unrotated)
            L_rotated <- vm$loadings # This is class 'loadings', behaves like matrix
            # Convert to matrix strictly
            L_rotated <- matrix(as.numeric(L_rotated), nrow = nrow(L_rotated), ncol = ncol(L_rotated))

            # Rotated Scores = X_scaled %*% L_rotated %*% solve(t(L_rotated) %*% L_rotated) ?
            # Or easier: Scores = X_scaled %*% solve(Correlation) %*% L_rotated (Regression method)
            # For approximate component scores:
            # Rot = vm$rotmat
            # Scores_rotated = Scores_unrotated %*% Rot
            S_unrotated <- U * scale_factor # = X_scaled %*% V (if standardized)
            S_rotated <- S_unrotated %*% vm$rotmat

            colnames(L_rotated) <- paste0("Factor", 1:factors)
            colnames(S_rotated) <- paste0("Factor", 1:factors)

            # 5. Communality and Specific Variance
            # h2 = sum of squared loadings
            h2 <- rowSums(L_rotated^2)
            u2 <- 1 - h2 # Unique variance / Specific variance

            # Store Results
            private$fa_results <- list(
                Loadings = L_rotated,
                Scores = S_rotated,
                Communality = h2,
                SpecificVar = u2,
                BlockIndices = final_indices,
                Rotation = vm$rotmat,
                Pos = (self$blocks$Start[final_indices] + self$blocks$End[final_indices]) / 2
            )

            message("Factor Analysis complete. Retained ", factors, " factors.")
        },
        #' @description
        #' Reconstruct model-implied correlation matrix from FA results.
        #' G_smooth = Loadings %*% t(Loadings) + diag(SpecificVar)
        #' @return Implied correlation matrix (N_blocks x N_blocks subset)
        get_model_correlation = function() {
            if (is.null(private$fa_results)) stop("Run analyze_block_structure() first.")

            L <- private$fa_results$Loadings
            u2 <- private$fa_results$SpecificVar

            G_smooth <- tcrossprod(L)
            diag(G_smooth) <- diag(G_smooth) + u2

            return(G_smooth)
        },
        #' @description
        #' Analyze structural drivers of factor loadings.
        #' Regresses Communality ~ Variance + Position.
        #' @return A summary of the linear model.
        analyze_structural_drivers = function() {
            if (is.null(private$fa_results)) stop("Run analyze_block_structure() first.")
            if (is.null(self$significance)) stop("Significance/Variance table not found.")

            idx <- private$fa_results$BlockIndices

            # Build dataframe
            df <- data.frame(
                Comm = private$fa_results$Communality,
                # Safe access to variances
                Var = self$significance$Variance[idx],
                # Position (already stored in fa_results)
                Pos = private$fa_results$Pos
            )

            message("Testing structural drivers (Communality ~ Variance + Position)...")
            m1 <- lm(Comm ~ Var + Pos, data = df)

            return(summary(m1))
        },
        #' @description
        #' Compute Haplotype Relationship Matrix.
        #' @param n_cores Number of cores to use.
        compute_hrm = function(n_cores = 1) {
            if (is.null(self$haplo_mat)) stop("Haplotypes not encoded.")

            message("Computing Haplotype Relationship Matrix (HRM) (Parallel: ", n_cores > 1, ")...")

            n_samples <- nrow(self$haplo_mat)
            n_blocks <- ncol(self$haplo_mat)

            if (n_cores > 1) {
                future::plan(future::multisession, workers = n_cores)
            } else {
                future::plan(future::sequential)
            }

            # Map-Reduce approach
            # Split blocks into chunks
            if (n_cores > 1) {
                chunk_indices <- split(1:n_blocks, cut(1:n_blocks, n_cores, labels = FALSE))
            } else {
                chunk_indices <- list(1:n_blocks)
            }

            # Map: Compute partial K for each chunk
            partial_Ks <- future.apply::future_lapply(chunk_indices, function(indices) {
                K_partial <- matrix(0, n_samples, n_samples)
                for (b in indices) {
                    ids <- self$haplo_mat[, b]
                    K_partial <- K_partial + outer(ids, ids, `==`)
                }
                return(K_partial)
            }, future.seed = TRUE)

            # Reduce: Sum partial Ks
            K <- Reduce(`+`, partial_Ks)

            # Normalize by number of blocks
            K <- K / n_blocks

            self$hrm <- K
            message("HRM computation complete.")
        },
        #' @description
        #' Fit Kernel Ridge Regression model.
        #' @param lambda Regularization parameter.
        #' @param use_cg Whether to use Conjugate Gradient solver.
        #' @param tol Tolerance for CG.
        #' @param max_iter Maximum iterations for CG.
        fit_krr = function(lambda = 0.1, use_cg = NULL, tol = 1e-5, max_iter = 1000) {
            if (is.null(self$hrm)) stop("HRM not computed.")
            if (is.null(self$pheno)) stop("Phenotypes not loaded.")

            K <- self$hrm
            y <- self$pheno
            n <- nrow(K)

            # Heuristic: Use CG if N > 5000 or explicitly requested
            if (is.null(use_cg)) {
                use_cg <- n > 5000
            }

            if (use_cg) {
                message("Fitting KRR using Conjugate Gradient (N=", n, ")...")

                # Solve (K + lambda*I) alpha = y
                # A = K + lambda*I
                # A is symmetric positive definite

                # Initial guess
                alpha <- rep(0, n)
                r <- y - (K %*% alpha + lambda * alpha)
                p <- r
                rsold <- sum(r * r)

                for (i in 1:max_iter) {
                    # Ap = (K + lambda*I) p
                    Ap <- K %*% p + lambda * p

                    alpha_step <- rsold / sum(p * Ap)
                    alpha <- alpha + alpha_step * p
                    r <- r - alpha_step * Ap

                    rsnew <- sum(r * r)
                    if (sqrt(rsnew) < tol) {
                        message("CG converged at iteration ", i)
                        break
                    }

                    p <- r + (rsnew / rsold) * p
                    rsold <- rsnew
                }

                if (i == max_iter) warning("CG did not converge within max_iter")
            } else {
                message("Fitting KRR using direct inversion (N=", n, ")...")
                I <- diag(n)
                alpha <- solve(K + lambda * I, y)
            }

            return(alpha)
        },
        #' @description
        #' Cross-validate KRR model.
        #' @param k Number of folds.
        #' @param lambdas Vector of lambdas to test.
        #' @param n_cores Number of cores to use.
        #' @param folds Optional vector of fold assignments.
        cross_validate = function(k = 5, lambdas = NULL, n_cores = 1, folds = NULL) {
            if (is.null(self$hrm)) stop("HRM not computed.")
            if (is.null(self$pheno)) stop("Phenotypes not loaded.")

            if (is.null(lambdas)) {
                lambdas <- 10^seq(-5, 5, length.out = 20)
            }

            n <- nrow(self$hrm)
            if (is.null(folds)) {
                folds <- sample(rep(1:k, length.out = n))
            } else {
                if (length(folds) != n) stop("Length of folds must match number of samples.")
                k <- length(unique(folds))
            }

            message("Running ", k, "-fold Cross-Validation on ", length(lambdas), " lambdas (Parallel: ", n_cores > 1, ")...")

            if (n_cores > 1) {
                future::plan(future::multisession, workers = n_cores)
            } else {
                future::plan(future::sequential)
            }

            # Grid search
            # Parallelize over lambdas? Or folds?
            # Parallelizing over lambdas is easier to aggregate.

            progressr::with_progress({
                p <- progressr::progressor(steps = length(lambdas))

                results <- future.apply::future_lapply(lambdas, function(lam) {
                    p()
                    mse_sum <- 0

                    for (i in 1:k) {
                        test_idx <- which(folds == i)
                        train_idx <- which(folds != i)

                        # Split K
                        K_train <- self$hrm[train_idx, train_idx]
                        K_test_train <- self$hrm[test_idx, train_idx]
                        y_train <- self$pheno[train_idx]
                        y_test <- self$pheno[test_idx]

                        # Train (using direct solve for stability in CV, usually N is smaller)
                        # Or use CG if N is large? Let's stick to solve for robustness in CV
                        # unless N_train is huge.

                        n_train <- length(train_idx)
                        I <- diag(n_train)

                        # Try-catch for singular matrices
                        tryCatch(
                            {
                                alpha <- solve(K_train + lam * I, y_train)
                                y_pred <- K_test_train %*% alpha
                                mse_sum <- mse_sum + sum((y_test - y_pred)^2)
                            },
                            error = function(e) {
                                mse_sum <- mse_sum + Inf
                            }
                        )
                    }

                    return(mse_sum / n)
                }, future.seed = TRUE)
            })

            mses <- unlist(results)
            best_idx <- which.min(mses)
            best_lambda <- lambdas[best_idx]

            message("CV Complete. Best lambda: ", best_lambda, " (MSE: ", mses[best_idx], ")")

            return(list(
                best_lambda = best_lambda,
                lambdas = lambdas,
                mse = mses
            ))
        },
        #' @description
        #' Plot Manhattan plot of local GEBV significance.
        #' @param threshold Significance threshold (p-value).
        plot_manhattan = function(threshold = 0.05) {
            if (is.null(self$significance)) stop("Significance testing not run.")
            if (is.null(self$blocks)) stop("Blocks not defined.")

            # Prepare data for plotting
            df <- data.frame(
                Block = self$blocks$BlockID,
                Start = self$blocks$Start,
                End = self$blocks$End,
                P_Value = self$significance$P_Value
            )

            # Calculate midpoint for plotting
            df$Pos <- (df$Start + df$End) / 2

            # Log p-values
            df$logP <- -log10(df$P_Value)

            # Threshold line
            thresh_log <- -log10(threshold)

            # Determine colors
            cols <- "black"
            if (!is.null(self$map) && "chr" %in% names(self$map)) {
                # Map blocks to chromosomes
                # We need to know which chromosome each block belongs to.
                # We can take the chromosome of the start marker of the block.

                # Ensure map is loaded and has chr
                real_start <- if (!is.null(self$active_markers)) self$active_markers[df$Start] else df$Start
                block_chrs <- self$map$chr[real_start]

                # Handle non-numeric chromosomes (e.g. "X", "Y") by converting to factor then integer
                if (!is.numeric(block_chrs)) {
                    block_chrs <- as.integer(as.factor(block_chrs))
                }

                # Alternating colors: black, grey
                # 1 -> black (1), 2 -> grey (2), 3 -> black (1)...
                # (chr %% 2) gives 1 or 0. +1 gives 2 or 1.
                col_vec <- c("grey", "black") # 0+1=1 -> grey, 1+1=2 -> black.
                # Usually odd=black, even=grey.
                # If chr=1, 1%%2 = 1. index 2.
                # If chr=2, 2%%2 = 0. index 1.

                cols <- col_vec[(block_chrs %% 2) + 1]
            }

            # Base R Plot
            plot(df$Pos, df$logP,
                pch = 20,
                col = cols,
                xlab = "Genomic Position (Index)",
                ylab = "-log10(P-value)",
                main = "Manhattan Plot (Local GEBV)",
                las = 1
            )

            abline(h = thresh_log, col = "red", lty = 2)

            # Return data invisibly instead of a plot object
            invisible(df)
        },
        #' @description
        #' Plot PCA of HRM.
        #' @param groups Optional vector of groups for coloring.
        plot_pca = function(groups = NULL) {
            if (is.null(self$hrm)) stop("HRM not computed. Run compute_hrm() first.")

            message("Calculating Eigen decomposition of HRM...")
            eig <- eigen(self$hrm, symmetric = TRUE)

            # Extract top 2 PCs
            pc1 <- eig$vectors[, 1]
            pc2 <- eig$vectors[, 2]

            # Variance explained
            var_expl <- eig$values / sum(eig$values) * 100

            # Plot
            xlab <- paste0("PC1 (", round(var_expl[1], 1), "%)")
            ylab <- paste0("PC2 (", round(var_expl[2], 1), "%)")

            col_vec <- "black"
            if (!is.null(groups)) {
                if (length(groups) != nrow(self$hrm)) {
                    warning("Length of groups does not match number of individuals. Ignoring groups.")
                } else {
                    # Convert groups to factor to get integer codes for colors
                    groups <- as.factor(groups)
                    col_vec <- as.integer(groups)
                }
            }

            plot(pc1, pc2,
                xlab = xlab, ylab = ylab,
                main = "PCA of Haplotype Relationship Matrix",
                pch = 19,
                col = col_vec
            )

            if (!is.null(groups)) {
                legend("topright", legend = levels(groups), col = 1:length(levels(groups)), pch = 19)
            }
        },
        #' @description
        #' Plot heatmap of local GEBVs.
        #' @param block_range Optional range of blocks to plot.
        plot_gebv_image = function(block_range = NULL) {
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated.")

            mat <- self$local_gebv$matrix

            if (!is.null(block_range)) {
                if (max(block_range) > ncol(mat)) stop("Block range exceeds number of blocks.")
                mat <- mat[, block_range, drop = FALSE]
            }

            # image() expects x as rows and y as cols.
            # We want X-axis = Blocks (Columns of mat), Y-axis = Individuals (Rows of mat).
            # image(x, y, z)
            # z should be matrix of dim length(x) * length(y).
            # If we pass mat directly: rows are x, cols are y.
            # So X-axis would be Individuals, Y-axis would be Blocks.
            # We want X-axis = Blocks. So we need to transpose.

            # t(mat) has dimensions (Blocks x Individuals).
            # Rows = Blocks, Cols = Individuals.
            # So X-axis = Blocks, Y-axis = Individuals.

            # Also, image() draws row 1 at bottom.
            # Usually we want Individual 1 at top? Or bottom is fine.

            message("Plotting GEBV heatmap...")

            # Palette
            pal <- colorRampPalette(c("blue", "white", "red"))(100)

            # Transpose for image
            # t(mat)

            image(
                x = 1:ncol(mat),
                y = 1:nrow(mat),
                z = t(mat),
                col = pal,
                xlab = "Block Index",
                ylab = "Individual Index",
                main = "Local GEBV Heatmap",
                useRaster = TRUE
            ) # Faster for large matrices
        },
        #' @description
        #' Identify superior haplotypes.
        #' @param top_n Number of top blocks to consider.
        identify_superior_haplotypes = function(top_n = 50) {
            if (is.null(self$significance)) stop("Significance testing not run.")
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated.")
            if (is.null(self$haplo_mat)) stop("Haplotypes not encoded.")

            message("Identifying superior haplotypes for top ", top_n, " blocks...")

            # Sort blocks by significance (p-value ascending)
            sorted_sig <- self$significance[order(self$significance$P_Value), ]
            top_blocks <- head(sorted_sig, top_n)

            results <- list()

            for (i in 1:nrow(top_blocks)) {
                block_id <- top_blocks$BlockID[i]

                # Get local GEBVs for this block
                # local_gebv$matrix is (Individuals x Blocks)
                gebvs <- self$local_gebv$matrix[, block_id]

                # Get haplotype IDs for this block
                haplos <- self$haplo_mat[, block_id]

                # Calculate mean GEBV per haplotype
                # Aggregate
                agg <- aggregate(gebvs, by = list(haplos), FUN = mean)
                colnames(agg) <- c("HaploID", "MeanGEBV")

                # Find best haplotype (highest positive effect)
                best_idx <- which.max(agg$MeanGEBV)
                best_haplo <- agg$HaploID[best_idx]
                effect_size <- agg$MeanGEBV[best_idx]

                results[[i]] <- data.table::data.table(
                    BlockID = block_id,
                    BestHaploID = best_haplo,
                    EffectSize = effect_size
                )
            }

            res_dt <- data.table::rbindlist(results)
            return(res_dt)
        },
        #' @description
        #' Calculate stacking scores.
        #' @param superior_haplos Data.table of superior haplotypes.
        score_stacking = function(superior_haplos) {
            if (is.null(self$haplo_mat)) stop("Haplotypes not encoded.")

            message("Scoring individuals based on superior haplotypes...")

            n_samples <- nrow(self$haplo_mat)
            scores <- rep(0, n_samples)

            # Iterate through superior haplotypes
            for (i in 1:nrow(superior_haplos)) {
                block_id <- superior_haplos$BlockID[i]
                best_id <- superior_haplos$BestHaploID[i]

                # Check which individuals have this haplotype
                # haplo_mat is (Individuals x Blocks)
                matches <- self$haplo_mat[, block_id] == best_id

                scores <- scores + as.integer(matches)
            }

            return(scores)
        },
        #' @description
        #' Plot stacking trend.
        #' @param scores Vector of stacking scores.
        #' @param superior_haplos Optional table of superior haplotypes (if scores not provided).
        plot_stacking_trend = function(scores = NULL, superior_haplos = NULL) {
            if (is.null(self$pheno)) stop("Phenotypes not loaded.")

            # If scores not provided, try to compute them
            if (is.null(scores)) {
                if (is.null(superior_haplos)) {
                    # Try to find top haplotypes automatically?
                    # Let's default to top 50 if nothing provided
                    message("No scores or superior haplotypes provided. Identifying top 50...")
                    superior_haplos <- self$identify_superior_haplotypes(top_n = 50)
                }
                scores <- self$score_stacking(superior_haplos)
            }

            message("Plotting stacking trend...")

            # Plot
            plot(scores, self$pheno,
                xlab = "Stacking Score (Number of Superior Haplotypes)",
                ylab = "Phenotype",
                main = "Stacking Validation",
                pch = 19,
                col = rgb(0, 0, 0, 0.5)
            ) # Semi-transparent black

            # Add regression line
            fit <- lm(self$pheno ~ scores)
            abline(fit, col = "red", lwd = 2)

            # Add R2 to plot
            r2 <- summary(fit)$r.squared
            legend("topleft", legend = paste0("R2 = ", round(r2, 3)), bty = "n")

            invisible(scores)
        },

        #' @description
        #' Impute missing genotypes.
        #' @param method Imputation method. Currently only "mean" is supported.
        impute_genotypes = function(method = "mean") {
            if (is.null(self$geno)) stop("Genotypes not loaded.")

            message("Imputing missing genotypes (method=", method, ")...")

            if (method == "mean") {
                n_markers <- ncol(self$geno)
                block_size <- 1000

                for (i in seq(1, n_markers, by = block_size)) {
                    ind <- i:min(i + block_size - 1, n_markers)

                    # Read chunk
                    chunk <- self$geno[, ind]

                    # Impute in memory
                    if (any(is.na(chunk))) {
                        # Column means
                        mus <- colMeans(chunk, na.rm = TRUE)

                        for (j in 1:ncol(chunk)) {
                            if (any(is.na(chunk[, j]))) {
                                # Round to nearest integer (0, 1, 2)
                                chunk[is.na(chunk[, j]), j] <- round(mus[j])
                            }
                        }

                        # Write back
                        self$geno[, ind] <- chunk
                    }
                }
            } else {
                stop("Only 'mean' imputation is currently supported for double FBMs.")
            }

            message("Imputation complete.")
        },

        #' @description
        #' Filter monomorphic markers (zero variance).
        filter_monomorphic = function() {
            if (is.null(self$geno)) stop("Genotypes not loaded.")

            message("Filtering monomorphic markers...")

            # Calculate standard deviations
            stats <- bigstatsr::big_scale()(self$geno)
            sds <- stats$scale

            # Identify non-zero variance markers
            keep_idx <- which(sds > 1e-8)

            n_total <- ncol(self$geno)
            n_keep <- length(keep_idx)
            n_drop <- n_total - n_keep

            if (n_drop > 0) {
                message("Dropped ", n_drop, " monomorphic markers.")
                self$active_markers <- keep_idx
            } else {
                message("No monomorphic markers found.")
                self$active_markers <- 1:n_total
            }
        },

        #' @description
        #' Save the project to an RDS file.
        #' @param path Path to the output .rds file.
        save_project = function(path) {
            if (is.null(self$geno)) stop("Nothing to save.")

            # Ensure path ends with .rds
            if (!endsWith(path, ".rds")) path <- paste0(path, ".rds")

            message("Saving project to ", path, "...")

            # Save the R6 object
            saveRDS(self, path)

            # Check if backing file is in the same directory
            bk_file <- paste0(private$backing_file, ".bk")
            if (file.exists(bk_file)) {
                message("Note: The backing file (.bk) is located at: ", bk_file)
                message("Ensure this file is kept with the .rds file if moving the project.")
            }
        },

        #' @description
        #' Scale Haplotype Effects.
        #' Centers the local GEBV matrix columns to have a mean of zero.
        #' @param local_gebv A matrix of local GEBVs (N x n_blocks). If NULL, uses self$local_gebv.
        #' @return A matrix of scaled local GEBVs.
        scale_haplo_effects = function(local_gebv = NULL) {
            if (is.null(local_gebv)) {
                if (is.null(self$local_gebv)) stop("Local GEBVs not loaded or computed.")
                local_gebv <- self$local_gebv$matrix
            } else if (is.list(local_gebv) && !is.data.frame(local_gebv) && "matrix" %in% names(local_gebv)) {
                local_gebv <- local_gebv$matrix
            }
            # Center columns to mean 0
            scaled_gebv <- scale(local_gebv, center = TRUE, scale = FALSE)
            return(scaled_gebv)
        },

        #' @description
        #' Analyze Haplotypes of Interest (HOI).
        #' Identifies superior haplotypes at a specific block by analyzing distribution.
        #' @param block_id The ID of the block to analyze.
        #' @param local_gebv A matrix of local GEBVs. If NULL, uses self$local_gebv.
        #' @return A list containing peak values, p-value, HOI haplotypes, and stats.
        analyze_hoi = function(block_id, local_gebv = NULL) {
            if (is.null(self$haplo_mat)) stop("Haplotypes must be encoded first.")
            if (is.null(local_gebv)) {
                if (is.null(self$local_gebv)) stop("Local GEBVs not loaded or computed.")
                local_gebv <- self$local_gebv$matrix
            } else if (is.list(local_gebv) && !is.data.frame(local_gebv) && "matrix" %in% names(local_gebv)) {
                local_gebv <- local_gebv$matrix
            }

            # Check if block_id is valid
            if (!block_id %in% 1:ncol(local_gebv)) stop("Invalid block_id.")

            effects <- local_gebv[, block_id]

            # Scale effects using internal method (or direct scale)
            # scaled_effects <- self$scale_haplo_effects(matrix(effects, ncol=1)) # Returns matrix
            # Use direct scale to return vector/single column matrix without complications
            scaled_effects <- scale(effects, center = TRUE, scale = FALSE)

            # Estimate density
            d <- stats::density(scaled_effects)

            # Find peaks (local maxima)
            peak_indices <- which(diff(sign(diff(d$y))) == -2) + 1
            peak_x <- d$x[peak_indices]
            peak_y <- d$y[peak_indices]

            if (length(peak_indices) < 2) {
                warning("Less than 2 peaks detected. Returning simple stats.")
                return(list(
                    peaks = peak_x,
                    p_value = NA,
                    nadir = mean(scaled_effects, na.rm = TRUE),
                    hoi_haplotypes = NULL,
                    stats = summary(scaled_effects)
                ))
            }

            # Sort peaks by x value
            sorted_peaks <- sort(peak_x)

            # Identify the two most prominent peaks
            top_peaks_idx <- order(peak_y, decreasing = TRUE)[1:2]
            main_peaks_x <- sort(peak_x[top_peaks_idx])

            low_peak <- main_peaks_x[1]
            high_peak <- main_peaks_x[2]

            # Find nadir between them
            range_indices <- which(d$x > low_peak & d$x < high_peak)
            if (length(range_indices) > 0) {
                nadir_idx <- range_indices[which.min(d$y[range_indices])]
                nadir_x <- d$x[nadir_idx]
            } else {
                nadir_x <- (low_peak + high_peak) / 2
            }

            # Split individuals
            group_high <- which(scaled_effects > nadir_x)
            group_low <- which(scaled_effects <= nadir_x)

            # Significance test
            if (length(group_high) > 1 && length(group_low) > 1) {
                test_res <- stats::t.test(scaled_effects[group_high], scaled_effects[group_low], alternative = "greater")
                p_val <- test_res$p.value
            } else {
                p_val <- NA
            }

            # Identify haplotypes associated with the high group
            haplo_col <- self$haplo_mat[, block_id]
            haplo_means <- tapply(scaled_effects, haplo_col, mean)
            hoi_haplos <- as.integer(names(haplo_means[haplo_means > nadir_x]))

            return(list(
                peaks = main_peaks_x,
                nadir = nadir_x,
                p_value = p_val,
                hoi_haplotypes = hoi_haplos,
                stats = list(
                    mean_high = mean(scaled_effects[group_high]),
                    mean_low = mean(scaled_effects[group_low]),
                    n_high = length(group_high),
                    n_low = length(group_low)
                )
            ))
        },

        #' @description
        #' Generate a Genomic Architecture Biplot.
        #' Plots individuals based on HRM PCA and overlays vectors for high-variance haploblocks.
        #' @param top_n Number of high-variance blocks to display as vectors.
        #' @param groups Optional vector of groups for coloring individuals.
        #' @param scale_vectors Scalar to adjust the length of vectors for visualization. If NULL, auto-scales.
        #' @param label_blocks Boolean, whether to label the block vectors with their IDs.
        #' @param highlight_ind Optional. Vector of sample names or indices to highlight.
        #' @description
        #' Generate a Genomic Architecture Biplot (Base R).
        #' Plots individuals based on HRM PCA and overlays vectors for high-variance haploblocks.
        #' @param top_n Number of high-variance blocks to display as vectors.
        #' @param groups Optional vector of groups for coloring individuals.
        #' @param scale_vectors Scalar to adjust the length of vectors for visualization. If NULL, auto-scales.
        #' @param label_blocks Boolean, whether to label the block vectors with their IDs.
        #' @param highlight_ind Optional. Vector of sample names or indices to highlight.
        plot_haplo_biplot = function(top_n = 10, groups = NULL, scale_vectors = NULL, label_blocks = TRUE, highlight_ind = NULL) {
            # 1. Validation
            if (is.null(self$hrm)) stop("HRM not computed. Run compute_hrm() first.")
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated. Run calculate_local_gebv() first.")
            if (is.null(self$significance)) stop("Significance/Variance not calculated. Run test_significance() first.")

            # 2. Get PC Scores (Individuals)
            message("Calculating PCA of Centered HRM...")

            # Helper to get centered HRM
            K <- self$hrm
            n <- nrow(K)
            row_means <- rowMeans(K)
            mean_K <- mean(K)

            # Use sweep for broadcasting subtract; K_centered = (I - 1/n J) K (I - 1/n J)
            Kc <- sweep(K, 1, row_means, "-")
            Kc <- sweep(Kc, 2, row_means, "-")
            Kc <- Kc + mean_K

            eig <- eigen(Kc, symmetric = TRUE)
            pc_scores <- data.frame(
                PC1 = eig$vectors[, 1],
                PC2 = eig$vectors[, 2]
            )

            # Calculate variance explained
            var_expl <- eig$values / sum(eig$values) * 100
            var_expl <- pmax(var_expl, 0) # Safety

            lab_x <- paste0("PC1 (", round(var_expl[1], 1), "%)")
            lab_y <- paste0("PC2 (", round(var_expl[2], 1), "%)")

            # 3. Loadings from Top Blocks
            # If fa_results exists, maybe use those? Protocol says "Top Variance Blocks" for biplot.

            # Sort blocks by Variance
            top_blocks_idx <- order(self$significance$Variance, decreasing = TRUE)[1:top_n]
            top_blocks_ids <- self$significance$BlockID[top_blocks_idx]

            local_gebv_mat <- if (is.list(self$local_gebv) && "matrix" %in% names(self$local_gebv)) self$local_gebv$matrix else self$local_gebv
            sel_gebvs <- local_gebv_mat[, top_blocks_idx, drop = FALSE]

            # Correlation with PC scores
            loadings <- data.frame(
                BlockID = top_blocks_ids,
                v1 = cor(sel_gebvs, pc_scores$PC1),
                v2 = cor(sel_gebvs, pc_scores$PC2)
            )

            # 4. Scaling
            max_pc <- max(abs(c(pc_scores$PC1, pc_scores$PC2)))
            max_load <- max(abs(c(loadings$v1, loadings$v2)))

            if (is.null(scale_vectors)) {
                if (is.na(max_load) || max_load == 0) max_load <- 1
                scaling_factor <- max_pc / max_load * 0.8
            } else {
                scaling_factor <- scale_vectors
            }

            loadings$v1 <- loadings$v1 * scaling_factor
            loadings$v2 <- loadings$v2 * scaling_factor

            # 5. Plotting
            # Setup colors
            if (is.null(groups)) {
                pt_col <- rgb(0, 0, 0, 0.3)
                pt_pch <- 16
            } else {
                if (length(groups) != n) {
                    warning("Length of groups does not match n_samples.")
                    pt_col <- rgb(0, 0, 0, 0.3)
                    pt_pch <- 16
                } else {
                    u_groups <- unique(groups)
                    # Simple color palette logic
                    pal <- c("red", "blue", "green", "orange", "purple", "cyan", "brown")
                    if (length(u_groups) > length(pal)) pal <- rainbow(length(u_groups))
                    pt_col <- pal[match(groups, u_groups)]
                    pt_pch <- 16
                }
            }

            x_lim <- range(c(pc_scores$PC1, loadings$v1), na.rm = TRUE)
            y_lim <- range(c(pc_scores$PC2, loadings$v2), na.rm = TRUE)

            # Add margin
            x_lim <- x_lim + c(-1, 1) * diff(x_lim) * 0.1
            y_lim <- y_lim + c(-1, 1) * diff(y_lim) * 0.1

            plot(pc_scores$PC1, pc_scores$PC2,
                xlab = lab_x, ylab = lab_y,
                main = "Genomic Architecture Biplot",
                col = pt_col, pch = pt_pch,
                xlim = x_lim, ylim = y_lim
            )

            grid()
            abline(h = 0, v = 0, lty = 2, col = "grey")

            # Vectors
            arrows(0, 0, loadings$v1, loadings$v2, col = "darkred", lwd = 1.5, length = 0.1)

            if (label_blocks) {
                text(loadings$v1, loadings$v2,
                    labels = paste0("B", loadings$BlockID),
                    pos = 1, col = "darkred", cex = 0.8
                )
            }

            # Highlight checks
            if (!is.null(highlight_ind)) {
                # Simple index matching
                idx <- NULL
                lbl <- NULL

                if (is.numeric(highlight_ind)) {
                    idx <- highlight_ind
                    lbl <- paste0("C", idx)
                } else if (is.character(highlight_ind)) {
                    if (!is.null(self$sample_ids)) {
                        idx <- match(highlight_ind, self$sample_ids)
                        lbl <- highlight_ind
                    } else if (!is.null(rownames(self$hrm))) {
                        idx <- match(highlight_ind, rownames(self$hrm))
                        lbl <- highlight_ind
                    }
                }

                if (!is.null(idx)) {
                    valid <- !is.na(idx)
                    points(pc_scores$PC1[idx[valid]], pc_scores$PC2[idx[valid]], col = "blue", pch = 17, cex = 1.5)
                    text(pc_scores$PC1[idx[valid]], pc_scores$PC2[idx[valid]], labels = lbl[valid], pos = 3, col = "blue", font = 2)
                }
            }

            invisible(list(scores = pc_scores, loadings = loadings))
        },
        #' @description
        #' Plot Factor Analysis Results (Stacked Manhattan Plots).
        #' @param factors Vector of factor indices to plot (default: all).
        plot_fa_genome = function(factors = NULL) {
            if (is.null(private$fa_results)) stop("Run analyze_block_structure() first.")

            res <- private$fa_results
            L <- res$Loadings
            Pos <- res$Pos

            n_factors <- ncol(L)
            if (is.null(factors)) factors <- 1:n_factors

            # Setup layout
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))

            par(mfrow = c(length(factors), 1), mar = c(2, 4, 2, 1), oma = c(2, 0, 2, 0))

            for (i in factors) {
                if (i > n_factors) next

                load_vec <- L[, i]
                # Color positive/negative
                cols <- ifelse(load_vec > 0, "blue", "red")

                plot(Pos, load_vec,
                    type = "h", lwd = 1.5, col = cols,
                    ylim = c(min(L), max(L)),
                    ylab = paste("Factor", i),
                    main = "",
                    xaxt = "n"
                ) # Suppress x-axis inside stack
                abline(h = 0)

                # Add simple axis if last
                if (i == tail(factors, 1)) {
                    axis(1)
                    mtext("Genomic Position", side = 1, line = 2.5, cex = 0.8)
                }
            }
            mtext("Latent Genomic Gradients", side = 3, outer = TRUE, line = 0.5, font = 2)
        },
        #' @description
        #' Plot Communality vs Position.
        #' Visualizes structural drivers (e.g. centromeres).
        plot_communality = function() {
            if (is.null(private$fa_results)) stop("Run analyze_block_structure() first.")

            comm <- private$fa_results$Communality
            pos <- private$fa_results$Pos

            plot(pos, comm,
                pch = 16, col = rgb(0, 0, 0, 0.4),
                xlab = "Genomic Position", ylab = "Communality (h2)",
                main = "Structural Constraint Map"
            )

            # Smooth trend if enough data and no NAs
            if (length(pos) > 5 && all(is.finite(pos)) && all(is.finite(comm))) {
                tryCatch(
                    {
                        lines(stats::loess.smooth(pos, comm, span = 0.5), col = "red", lwd = 2)
                    },
                    error = function(e) warning("Could not plot smooth trend: ", e$message)
                )
            }
            grid()
        },
        #' @description
        #' Plot Scree Plot of Factor Analysis.
        #' Shows % variance explaining by each factor or singular value.
        plot_scree = function() {
            if (is.null(private$fa_results)) stop("Run analyze_block_structure() first.")

            # We don't have all SVD values stored in fa_results currently,
            # only the rotated loadings for chosen factors.
            # Ideally analyze_block_structure should store eigenvalues from SVD.
            # But we can calculate variance explained by the extracted factors.

            L <- private$fa_results$Loadings
            # Sum of squared loadings per factor = Variance explained by that factor (after rotation)
            var_per_factor <- colSums(L^2)
            # Total variance = Number of items? (Since inputs were standardized)
            total_var <- nrow(L)

            pve <- var_per_factor / total_var * 100

            # Plot
            barplot(pve,
                names.arg = colnames(L),
                col = "steelblue",
                ylab = "% Variance Explained",
                main = "Scree Plot (Extracted Factors)",
                ylim = c(0, max(pve) * 1.2)
            )

            text(x = seq_along(pve) * 1.2 - 0.5, y = pve, labels = paste0(round(pve, 1), "%"), pos = 3, cex = 0.8)
        },

        #' @description
        #' Plot HOI Effect Distribution with Checks.
        #' Visualizes the distribution of scaled effects for a block and key check genotypes.
        #' @param block_id Integer. The ID of the block to visualize.
        #' @param highlight_ind Optional. Vector of sample names or indices to highlight as checks.
        plot_hoi_distribution = function(block_id, highlight_ind = NULL) {
            # 1. Validation
            if (is.null(self$local_gebv)) stop("Local GEBVs not calculated. Run calculate_local_gebv() first.")

            # 2. Get Data
            # local_gebv is list with $matrix
            local_gebv_mat <- if (is.list(self$local_gebv) && "matrix" %in% names(self$local_gebv)) self$local_gebv$matrix else self$local_gebv

            if (block_id < 1 || block_id > ncol(local_gebv_mat)) stop("Invalid block_id.")

            effects <- local_gebv_mat[, block_id]
            scaled_effects <- scale(effects, center = TRUE, scale = FALSE)

            # Run analysis to get nadir/peaks for context
            hoi_res <- self$analyze_hoi(block_id)
            nadir <- hoi_res$nadir

            # 3. Prepare Plot Data
            d <- stats::density(as.vector(scaled_effects))

            # Base R Plot
            plot(d,
                main = paste("Effect Distribution - Block", block_id),
                xlab = "Scaled Effect", ylab = "Density", col = "blue", lwd = 2
            )

            # Fill
            polygon(d, col = rgb(0.68, 0.85, 0.9, 0.5), border = "blue") # lightblue

            # Nadir
            grid()
            abline(v = nadir, col = "darkgreen", lty = 2, lwd = 2)
            mtext("Separation", at = nadir, col = "darkgreen", line = 0.5, cex = 0.8)

            # 4. Highlight Checks
            if (!is.null(highlight_ind)) {
                if (is.character(highlight_ind)) {
                    warning("Character names for checks might not be supported if matrix has no rownames. Using indices is safer.")
                }

                check_vals <- NULL
                check_labels <- NULL

                # Check for valid indices first
                if (is.numeric(highlight_ind)) {
                    check_vals <- scaled_effects[highlight_ind]
                    check_labels <- if (!is.null(names(highlight_ind))) names(highlight_ind) else paste0("Check:", highlight_ind)
                } else if (is.character(highlight_ind)) {
                    # Try to match against rownames if they exist
                    rn <- rownames(local_gebv_mat)
                    if (!is.null(rn)) {
                        match_idx <- match(highlight_ind, rn)
                        valid_match <- !is.na(match_idx)
                        if (any(valid_match)) {
                            check_vals <- scaled_effects[match_idx[valid_match]]
                            check_labels <- highlight_ind[valid_match]
                        }
                    }
                    # If names didn't work and we are desperate?
                    if (is.null(check_vals)) {
                        warning("Could not match character names to data rownames. Ensure matrix has rownames.")
                    }
                }

                # Filter valid
                if (!is.null(check_vals)) {
                    valid_idx <- which(is.finite(check_vals))
                    if (length(valid_idx) > 0) {
                        check_vals <- check_vals[valid_idx]
                        check_labels <- check_labels[valid_idx]

                        abline(v = check_vals, col = "red", lty = 1, lwd = 1.5)
                        # Add labels on top axis or above
                        axis(3, at = check_vals, labels = check_labels, col.axis = "red", col = "red", las = 2, cex.axis = 0.7)
                    } else {
                        warning("No finite values found for checks.")
                    }
                }
            }

            invisible(hoi_res)
        }
    ),
    private = list(
        #' @field backing_file Path to backing file.
        backing_file = NULL,
        #' @field fa_results List to store Factor Analysis results.
        fa_results = NULL,

        #' @description
        #' Read CSV/TXT genotypes with robust parsing.
        #' @param path Path to file.
        read_csv_genotypes = function(path) {
            message("Reading genotypes from ", path, "...")

            # Fast read
            dt <- data.table::fread(path, data.table = FALSE)

            # Heuristic: Check if first column is ID
            # If first column is character and unique, treat as ID and remove
            first_col <- dt[[1]]
            if (is.character(first_col) || is.factor(first_col)) {
                # Check uniqueness
                if (anyDuplicated(first_col) == 0) {
                    message("Detected Sample IDs in first column. Capturing as sample_ids.")
                    self$sample_ids <- first_col
                    dt <- dt[, -1, drop = FALSE]
                }
            }

            # Check for non-numeric data
            is_numeric <- all(vapply(dt, is.numeric, logical(1)))

            if (!is_numeric) {
                message("Detected non-numeric genotypes. Attempting to parse...")
                # Conversion function
                parse_geno <- function(x) {
                    if (is.numeric(x)) {
                        return(x)
                    }
                    x <- as.character(x)
                    # Replace common patterns
                    # 0/0, 0|0 -> 0
                    # 0/1, 0|1, 1/0, 1|0 -> 1
                    # 1/1, 1|1 -> 2
                    # . -> NA
                    res <- rep(NA_real_, length(x))

                    # Exact matches for speed
                    res[x %in% c("0", "0/0", "0|0", "A/A", "A|A")] <- 0
                    res[x %in% c("1", "0/1", "0|1", "1/0", "1|0", "A/B", "A|B", "B/A", "B|A")] <- 1
                    res[x %in% c("2", "1/1", "1|1", "B/B", "B|B")] <- 2

                    return(res)
                }

                # Apply to all columns
                # data.table set function is efficient but we converted to DF above.
                # Let's convert back to DT for fast processing if large?
                # Or just iterate since we're in-memory anyway before FBM.

                for (j in seq_len(ncol(dt))) {
                    if (!is.numeric(dt[[j]])) {
                        dt[[j]] <- parse_geno(dt[[j]])
                    }
                }
            }

            # Convert to matrix
            mat <- as.matrix(dt)

            # Check for NAs and warn
            n_na <- sum(is.na(mat))
            if (n_na > 0) {
                warning("Imported matrix contains ", n_na, " missing values (NA). Imputation may be needed.")
            }

            message("Converting to FBM...")
            self$geno <- bigstatsr::as_FBM(mat,
                type = "double",
                backingfile = private$backing_file
            )

            # Generate dummy map if needed to strictly satisfy bigsnpr-like expectations (though we use our own map usually)
            # If headers existed, they are colnames(dt)
            marker_names <- colnames(dt)
            if (length(marker_names) == ncol(self$geno)) {
                # Basic map
                self$map <- data.table::data.table(
                    marker = marker_names,
                    chr = 1,
                    pos = 1:ncol(self$geno),
                    ref = "A",
                    alt = "B"
                )
            }
        }
    )
)

#' Load HaploGeno Project
#'
#' @param path Path to the .rds file.
#' @return A HaploObject instance.
#' @export
load_haplo_project <- function(path) {
    if (!file.exists(path)) stop("File not found: ", path)

    message("Loading project from ", path, "...")

    # Load the R6 object
    project <- readRDS(path)

    # Re-attach/Reconstruct FBM
    if (!is.null(project$geno)) {
        # FBM objects lose their external pointer when saved/loaded via saveRDS
        # We need to reconstruct the FBM object to reconnect to the backing file.

        message("Re-connecting FBM backing file...")

        # Original backing file path
        original_bk <- project$geno$backingfile
        # Strip .bk extension if present
        original_bk_no_ext <- tools::file_path_sans_ext(original_bk)

        bk_path <- NULL

        # Check if .bk exists at original location
        if (file.exists(paste0(original_bk_no_ext, ".bk"))) {
            bk_path <- original_bk_no_ext
        } else {
            # Try same directory as the .rds file
            dir_name <- dirname(path)
            base_name <- basename(original_bk_no_ext)
            new_bk <- file.path(dir_name, base_name)

            if (file.exists(paste0(new_bk, ".bk"))) {
                bk_path <- new_bk
                message("Found backing file at new location: ", bk_path)
            } else {
                warning("Backing file not found. Genotypes will be inaccessible.")
            }
        }

        if (!is.null(bk_path)) {
            # Map integer type to string for FBM constructor
            type_map <- c(
                "1" = "unsigned char",
                "2" = "unsigned short",
                "4" = "integer",
                "6" = "float",
                "8" = "double"
            )
            type_code <- as.character(project$geno$type)
            if (type_code %in% names(type_map)) {
                type_str <- type_map[type_code]

                # Reconstruct FBM
                project$geno <- bigstatsr::FBM(
                    nrow = project$geno$nrow,
                    ncol = project$geno$ncol,
                    type = type_str,
                    backingfile = bk_path,
                    create_bk = FALSE,
                    is_read_only = project$geno$is_read_only
                )
            } else {
                warning("Unknown FBM type code: ", type_code, ". Cannot reconstruct FBM.")
            }
        }
    }

    return(project)
}
