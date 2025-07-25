#' Sample Causal Genetic Variants
#'
#' @description
#' Samples a specified number of causal variants from a set of GDS files. The function operates in one of three modes: uniform sampling, or weighted sampling based on LDAK or custom user-provided weights. It can filter sites by missingness or minor allele count (MAC) before sampling and returns the genotypes and annotations of the selected variants.
#'
#' @details
#' The sampling mode is determined automatically based on the provided arguments:
#' \itemize{
#'   \item \strong{uniform}: Default mode. Triggered if neither `ldak_ld_weights_paths` nor `sampling_weights_paths` is provided. Variants are sampled uniformly at random from the pool of available sites after filtering.
#'   \item \strong{ldak_weights}: Triggered by providing paths in `ldak_ld_weights_paths`. Requires chromosome information to be provided either within the files or via the `ldak_chrs` parameter.
#'   \item \strong{custom_weights}: Triggered by providing paths in `sampling_weights_paths`. Variants are sampled according to user-defined weights from files containing mandatory `chr`, `pos`, and `weight` columns.
#' }
#' The function robustly handles cases where the number of GDS files and weight files do not match, under the condition that one of the inputs is a single, consolidated file (e.g., one GDS file for all chromosomes, or one weight file for all chromosomes).
#'
#' The join between weights/annotations and genotypes is performed using a robust composite key (`chr`, `pos`, `ref`, `alt`) to ensure correctness.
#'
#' @param gds_paths A character vector of paths to the GDS files. When `length > 1`, each file is expected to contain data for a single chromosome.
#' @param n_causal_sites A single positive integer specifying the total number of causal sites to sample across all GDS files.
#' @param ldak_ld_weights_paths A character vector of paths to LDAK weights files. The file must have either two 
#'  or three columns. 
#' - Two-columns files: Raw (no header) *.short output file by `ldak --calc-weights`. The columns are `rsid` and `weight`.  
#' - Three-columns files: File with a header, containing columns `chr` `rsid` and `weight` (case-insensitive). 
#' See `ldak_chrs` for how chromosome information is handled for 2-column files. Cannot be used with `sampling_weights_paths`. 
#' @param sampling_weights_paths A character vector of paths to custom sampling weight files. Files must have a header and contain mandatory columns (case-insensitive): `chr`, `pos`, and `weight`. The columns `rsid`, `ref`, and `alt` are optional; a warning is issued if they are missing.
#' @param ldak_chrs A character vector used to assign chromosomes to the files specified in `ldak_ld_weights_paths`. This is required only when the LDAK files are in the standard 2-column format (rsid, weight) AND the number of LDAK weight files does not match the number of GDS files. Its length must match `length(ldak_ld_weights_paths)`. It is ignored if the LDAK files are in a 3-column format with a 'chr' column header.
#' @param stabilize_sampling_weights A logical flag. If `TRUE` (default), a small constant (10^-6) is added to all weights to prevent zero-probability sampling.
#' @param weights_power A numeric value used to raise the sampling weights to a power. Default is 1. For LDAK weights, consider `weights_power = -0.25` to enrich for variants in high-LD regions, or `weights_power = 0.25` for low-LD regions.
#' @param missingness_threshold A numeric value between 0 and 1. Variants with a missingness rate greater than this threshold will be excluded before sampling. Default is `NULL`.
#' @param mac_threshold A positive integer. Variants with a minor allele count (MAC) less than this threshold will be excluded before sampling. Default is `NULL`.
#' @param n_threads A positive integer specifying the number of threads for parallel operations on GDS files. Default is 1.
#'
#' @return
#' A list containing two data frames:
#' \describe{
#'   \item{`causal_annotation`}{A data frame with annotations for each of the `n_causal_sites` sampled variants.}
#'   \item{`causal_genotypes`}{A data frame with `n_causal_sites` rows and `(k + N)` columns, where `k` is the number of annotation columns and `N` is the number of samples. It contains the genotypes and annotations for the causal variants.}
#' }
#'
#' @importFrom data.table fread rbindlist fwrite
#' @importFrom dplyr %>% mutate select filter group_by summarise ungroup pull inner_join rename slice_sample cur_group
#' @importFrom glue glue
#' @importFrom gdsfmt showfile.gds
#' @importFrom SeqArray seqOpen seqGetData seqClose seqAlleleCount seqSetFilter seqMissing seqSummary
#' @importFrom SeqVarTools altDosage
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Ensure required packages and package data are available for the examples
#' if (!requireNamespace("data.table", quietly = TRUE) ||
#'     !requireNamespace("SeqArray", quietly = TRUE)) {
#'   stop("Please install data.table and SeqArray to run these examples.")
#' }
#'
#' ##================================================================
#' ## Example 1: Uniform sampling of causal sites
#' ##================================================================
#' cat(">>> Running Example 1: Uniform Sampling\n")
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#' n_causal_total <- 200
#' uniform_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_total,
#'   n_threads = 2
#' )
#' stopifnot(nrow(uniform_sites$causal_genotypes) == n_causal_total)
#' cat("OK: Uniform sampling returned the correct number of sites.\n")
#'
#' ##================================================================
#' ## Example 2: Weighted sampling with a custom weights file
#' ##================================================================
#' cat(">>> Running Example 2: Custom Weights File\n")
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#'
#' # Create temporary custom weight files
#' tmp_custom_files <- sapply(1:3, function(chr_idx) {
#'   tmp_file <- tempfile(fileext = ".tsv")
#'   gds <- SeqArray::seqOpen(gds_paths[chr_idx])
#'   on.exit(SeqArray::seqClose(gds))
#'   custom_weights_df <- data.frame(
#'     chr = SeqArray::seqGetData(gds, "chromosome"),
#'     pos = SeqArray::seqGetData(gds, "position"),
#'     rsid = SeqArray::seqGetData(gds, "annotation/id"),
#'     ref = SeqArray::seqGetData(gds, "$ref"),
#'     alt = SeqArray::seqGetData(gds, "$alt"),
#'     weight = runif(SeqArray::seqSummary(gds)$num.variant, 1, 100)
#'   )
#'   data.table::fwrite(custom_weights_df, file = tmp_file, sep = "\t")
#'   return(tmp_file)
#' })
#'
#' n_causal_custom <- 50
#' custom_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_custom,
#'   sampling_weights_paths = tmp_custom_files
#' )
#'
#' stopifnot(nrow(custom_sites$causal_genotypes) == n_causal_custom)
#' cat("OK: Custom weighted sampling returned the correct number of sites.\n")
#' file.remove(tmp_custom_files)
#'
#' ##================================================================
#' ## Example 3: Weighted sampling using LDAK weights and the `ldak_chrs` parameter
#' ##================================================================
#' cat(">>> Running Example 3: LDAK Weights from R Object")
#' data(phenocause.ldak_weights)
#' temp_dir <- tempdir()
#' # Create standard 2-column LDAK weight files (no header)
#' temp_ldak_paths <- sapply(names(phenocause.ldak_weights), function(chr_name) {
#'   file_path <- file.path(temp_dir, paste0(chr_name, ".weights.short"))
#'   data.table::fwrite(
#'     phenocause.ldak_weights[[chr_name]],
#'     file = file_path, sep = " ", col.names = FALSE
#'   )
#'   return(file_path)
#' })
#'
#' n_causal_ldak <- 300
#' # Since the files have no 'chr' column, we must supply it via `ldak_chrs`
#' ldak_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_ldak,
#'   ldak_ld_weights_paths = temp_ldak_paths,
#'   ldak_chrs = c("1", "2", "3"), # Must correspond to the files
#'   n_threads = 2
#' )
#' stopifnot(nrow(ldak_sites$causal_genotypes) == n_causal_ldak)
#' cat("OK: LDAK weighted sampling returned the correct number of sites.\n")
#' file.remove(temp_ldak_paths)
#' }
sample_causal_sites <- function(gds_paths,
                              n_causal_sites,
                              ldak_ld_weights_paths = NULL,
                              sampling_weights_paths = NULL,
                              ldak_chrs = NULL,
                              stabilize_sampling_weights = TRUE,
                              weights_power = 1,
                              missingness_threshold = NULL,
                              mac_threshold = NULL,
                              n_threads = 1) {

    # Helper function to safely open and close GDS files:
    safely_read_gds <- function(path, fun) {
        gds <- tryCatch(SeqArray::seqOpen(path), error = function(e) stop(glue::glue("!!! Failed to open GDS file: {path}. Error: {e$message}")))
        on.exit(SeqArray::seqClose(gds))
        fun(gds)
    }

    ##=====================================================================================================
    ### SECTION 1: INPUT VALIDATION
    ##=====================================================================================================
    if (!is.numeric(n_causal_sites) || (n_causal_sites %% 1) != 0 || n_causal_sites <= 0) stop("!!! `n_causal_sites` must be a single, positive integer.")
    if (!is.numeric(weights_power) || length(weights_power) != 1) stop("!!! `weights_power` must be a single numeric value.")
    if (!is.null(missingness_threshold) && (!is.numeric(missingness_threshold) || missingness_threshold < 0 || missingness_threshold > 1)) stop("!!! `missingness_threshold` must be a single numeric value between 0 and 1.")
    if (!is.null(mac_threshold) && (!is.numeric(mac_threshold) || (mac_threshold %% 1) != 0 || mac_threshold <= 0)) stop("!!! `mac_threshold` must be a single, positive integer.")
    if (!is.numeric(n_threads) || (n_threads %% 1) != 0 || n_threads <= 0) stop("!!! `n_threads` must be a single, positive integer.")

    if (!is.null(ldak_ld_weights_paths) && !is.null(sampling_weights_paths)) stop("!!! Cannot provide both `ldak_ld_weights_paths` and `sampling_weights_paths`.")
    if (any(duplicated(gds_paths))) stop("!!! Found duplicate paths in `gds_paths`.")

    sampling_mode <- "uniform"
    weights_paths <- NULL
    if (!is.null(sampling_weights_paths)) {
        sampling_mode <- "custom_weights"
        weights_paths <- sampling_weights_paths
    }
    if (!is.null(ldak_ld_weights_paths)) {
        sampling_mode <- "ldak_weights"
        weights_paths <- ldak_ld_weights_paths
    }
    
    if (sum(file.exists(gds_paths)) != length(gds_paths)) stop("!!! One or more GDS files are missing.\n")
    if (!is.null(weights_paths)) {
        if (sum(file.exists(weights_paths)) != length(weights_paths)) stop("!!! One or more weight files are missing.\n")
        if (any(duplicated(weights_paths))) stop("!!! Found duplicate paths in the weights vector.")
        
        n_gds <- length(gds_paths)
        n_weights <- length(weights_paths)
        if (n_gds != n_weights && n_gds != 1 && n_weights != 1) {
            stop(glue::glue("!!! Mismatch in file counts. GDS files ({n_gds}) and weight files ({n_weights}) do not match, and neither is of length 1."))
        }
    }
    
    gds_chrom_list <- lapply(gds_paths, safely_read_gds, fun = function(gds) unique(SeqArray::seqGetData(gds, "chromosome")))

    if (length(gds_paths) > 1 && any(sapply(gds_chrom_list, length) > 1)) {
        stop("!!! When multiple GDS files are provided, each file must contain data for only one chromosome.")
    } else if (length(gds_paths) == 1 && length(gds_chrom_list[[1]]) > 1) {
        warning("Single GDS file contains multiple chromosomes. Sampling will proceed from all of them.")
    }

    sample_ids <- safely_read_gds(gds_paths[1], fun = function(gds) SeqArray::seqGetData(gds, "sample.id"))
    if (length(gds_paths) > 1) {
        for (i in 2:length(gds_paths)) {
            current_ids <- safely_read_gds(gds_paths[i], fun = function(gds) SeqArray::seqGetData(gds, "sample.id"))
            if (!identical(sample_ids, current_ids)) {
                stop("!!! Sample IDs are not consistent across all provided GDS files.")
            }
        }
    }
    
    ##=====================================================================================================
    ### SECTION 2: AGGREGATE AND VALIDATE WEIGHTS
    ##=====================================================================================================
    weights_df <- NULL
    if (sampling_mode != "uniform") {
        # cat("DEBUG: Aggregating and validating weight files...\n")
        
        if (sampling_mode == "custom_weights") {
            res <- lapply(seq_along(weights_paths), function(i) {
                p <- weights_paths[i]
                df <- data.table::fread(p, header = TRUE)
                names(df) <- tolower(names(df))
                mand_cols <- c("chr", "pos", "weight")
                if (!all(mand_cols %in% names(df))) stop(glue::glue("!!! Custom weight file {p} is missing one or more mandatory columns: chr, pos, weight."))
                opt_cols <- c("rsid", "ref", "alt")
                missing_opts <- setdiff(opt_cols, names(df))
                list(data = df, missing = missing_opts, path = p)
            })
            weights_list <- lapply(res, `[[`, "data")
            missing_info <- lapply(res, function(r) if(length(r$missing)>0) r else NULL)
            missing_info <- missing_info[!sapply(missing_info, is.null)]
            if (length(missing_info) > 0) {
                msg <- paste(sapply(missing_info, function(m) {
                    glue::glue("File '{basename(m$path)}' missing: {paste(m$missing, collapse=', ')}")
                }), collapse="\n")
                warning(glue::glue("Optional columns were missing from custom weight files:\n{msg}"))
            }
            weights_df <- data.table::rbindlist(weights_list, use.names = TRUE, fill = TRUE)

        # } else if (sampling_mode == "ldak_weights") {
        #     if (!is.null(ldak_chrs) && length(ldak_chrs) != length(ldak_ld_weights_paths)) stop("!!! `ldak_chrs` length must match `ldak_ld_weights_paths` length.")
            
        #     weights_list <- lapply(seq_along(ldak_ld_weights_paths), function(i) {
        #         p <- ldak_ld_weights_paths[i]
        #         header_df <- data.table::fread(p, nrows = 1, header=FALSE)
                
        #         if (ncol(header_df) == 3 && all(sapply(header_df, is.character))) { # Likely 3-col with header
        #             df <- data.table::fread(p, header = TRUE)
        #             names(df) <- tolower(names(df))
        #             if (!all(c("chr", "rsid", "weight") %in% names(df))) stop(glue::glue("!!! 3-column LDAK file {p} has incorrect headers. Expected (case-insensitive): 'chr', 'rsid', 'weight'."))
        #         } else { # Assume 2-col, no header
        #             if (is.null(ldak_chrs)) stop("!!! 2-column LDAK files provided, but `ldak_chrs` parameter is NULL. Chromosome info is missing.")
        #             df <- data.table::fread(p, header = FALSE)
        #             names(df) <- c("rsid", "weight")
        #             df$chr <- ldak_chrs[i]
        #         }
        #         return(df)
        #     })
        #     weights_df <- data.table::rbindlist(weights_list, use.names = TRUE, fill = TRUE)
        # }
        } else if (sampling_mode == "ldak_weights") {
            if (!is.null(ldak_chrs) && length(ldak_chrs) != length(ldak_ld_weights_paths)) stop("!!! `ldak_chrs` length must match `ldak_ld_weights_paths` length.")
            
            weights_list <- lapply(seq_along(ldak_ld_weights_paths), function(i) {
                p <- ldak_ld_weights_paths[i]
                header_df <- data.table::fread(p, nrows = 1, header=FALSE)
                
                if (ncol(header_df) == 3 && all(sapply(header_df, is.character))) { # Likely 3-col with header
                    df <- data.table::fread(p, header = TRUE)
                    names(df) <- tolower(names(df))
                    if (!all(c("chr", "rsid", "weight") %in% names(df))) stop(glue::glue("!!! 3-column LDAK file {p} has incorrect headers. Expected (case-insensitive): 'chr', 'rsid', 'weight'."))
                
                } else { # Assume 2-col, no header
                    if (length(ldak_ld_weights_paths) > 1 && 
                        length(gds_paths) > 1 && 
                        length(ldak_ld_weights_paths) == length(gds_paths)) {
                        # Infer the chromosome from the corresponding GDS file.
                        df <- data.table::fread(p, header = FALSE)
                        names(df) <- c("rsid", "weight")
                        df$chr <- gds_chrom_list[[i]] # pre-validated GDS list
                    } else {
                        # require ldak_chrs for all other 2-column cases.
                        if (is.null(ldak_chrs)) stop("!!! 2-column LDAK files provided, but chromosome info cannot be inferred. Use the `ldak_chrs` parameter.")
                        df <- data.table::fread(p, header = FALSE)
                        names(df) <- c("rsid", "weight")
                        df$chr <- ldak_chrs[i]
                    }
                }
                return(df)
            })
            weights_df <- data.table::rbindlist(weights_list, use.names = TRUE, fill = TRUE)
        }

        initial_rows <- nrow(weights_df)
        essential_cols <- intersect(c("chr", "pos", "ref", "alt", "rsid", "weight"), names(weights_df))
        weights_df <- na.omit(weights_df, cols = essential_cols)
        rows_removed <- initial_rows - nrow(weights_df)
        if (rows_removed > 0) warning(glue::glue("{rows_removed} rows with NA in essential columns were removed from the combined weights data."))
        
        if (!is.numeric(weights_df$weight)) stop("!!! 'weight' column must be numeric.")
        if (any(weights_df$weight <= 0)) stop("!!! All weights must be strictly positive (> 0).")
        
        weights_df$chr <- as.character(weights_df$chr)
        if ("pos" %in% names(weights_df)) weights_df$pos <- as.integer(weights_df$pos) 
        if ("rsid" %in% names(weights_df)) weights_df$rsid <- as.character(weights_df$rsid)
    }
    
    ##=====================================================================================================
    ### SECTION 3: DEFINE SAMPLING POOL AND SITES PER CHROMOSOME
    ##=====================================================================================================
    # cat("DEBUG: Defining sampling pool and sites per chromosome...\n")
    
    ## Build the master annotation table 
    variant_list <- lapply(gds_paths, function(p) {
        safely_read_gds(p, function(gds) {
        
            ## total variants in file
            n_var <- SeqArray::seqSummary(gds, verbose = FALSE)$num.variant
            if (n_var == 0) return(NULL)
    
            ## apply MAC / missingness filters
            keep <- rep(TRUE, n_var)
            if (!is.null(mac_threshold)) {
                cat(">>> Filtering by MAC...\n")
                keep <- keep & (SeqArray::seqAlleleCount(gds, minor = TRUE) >= mac_threshold)
            }
            if (!is.null(missingness_threshold)) {
                cat(">>> Filtering by missingness...\n")
                keep <- keep & (SeqArray::seqMissing(gds, per.variant = TRUE) <= missingness_threshold)
            }
            if (!any(keep)) return(NULL)
    
            SeqArray::seqSetFilter(gds, variant.sel = keep)
            data.frame(
                chr = as.character(SeqArray::seqGetData(gds, "chromosome")),
                pos = as.integer(SeqArray::seqGetData(gds, "position")),
                rsid = as.character(SeqArray::seqGetData(gds, "annotation/id")),
                ref = as.character(SeqArray::seqGetData(gds, "$ref")),
                alt = as.character(SeqArray::seqGetData(gds, "$alt")),
                variant.id.internal = SeqArray::seqGetData(gds, "variant.id"),
                gds_path = p,
                stringsAsFactors = FALSE
            )
        })
    })
    
    ## drop chromosomes with no qualifying variants
    variant_list <- Filter(Negate(is.null), variant_list)
    if (length(variant_list) == 0) {
        stop("!!! No variants available for sampling after applying MAC/missingness filters.")
    }
    
    ## rbindlist is robust to differing column order
    all_variants_annot <- data.table::rbindlist(variant_list, use.names = TRUE, fill = TRUE)
    
    if (is.null(all_variants_annot) || nrow(all_variants_annot) == 0) stop("!!! No variants available for sampling after applying filters.")

    if (sampling_mode != "uniform") {
        join_by_cols <- intersect(c("chr", "pos", "rsid", "ref", "alt"), names(weights_df))
        all_variants_annot <- dplyr::inner_join(all_variants_annot, weights_df, by = join_by_cols)
    } else {
        all_variants_annot$weight <- 1
    }
    
    if (nrow(all_variants_annot) < n_causal_sites) {
        stop(glue::glue("!!! Cannot sample {n_causal_sites} sites. Only {nrow(all_variants_annot)} variants are available after filtering and merging with weights."))
    }
    
    constant <- ifelse(stabilize_sampling_weights, 1e-6, 0)
    all_variants_annot$final_weight <- (all_variants_annot$weight + constant)^weights_power
    
    chrom_summary <- all_variants_annot %>%
        dplyr::group_by(chr) %>%
        dplyr::summarise(total_weight = sum(final_weight), available_sites = dplyr::n(), .groups = 'drop')
        
    chrom_summary$n_causal <- round(n_causal_sites * (chrom_summary$total_weight / sum(chrom_summary$total_weight)))
    
    # Correction Pass 1: Cap at available sites and redistribute excess
    excess <- sum(pmax(0, chrom_summary$n_causal - chrom_summary$available_sites))
    chrom_summary$n_causal <- pmin(chrom_summary$n_causal, chrom_summary$available_sites)
    if (excess > 0) {
        warning("! After the filtering steps, the number of candidate causal variants is lower than the requested number
                    of causal variants for at least one chromosome. 
                    The number of sampled causal variants will be re-distributed across chromosomes\n")
        eligible_chrs <- chrom_summary$n_causal < chrom_summary$available_sites
        if(any(eligible_chrs)) {
            prop_weights <- chrom_summary$total_weight[eligible_chrs]
            adj <- round(excess * (prop_weights / sum(prop_weights)))
            chrom_summary$n_causal[eligible_chrs] <- chrom_summary$n_causal[eligible_chrs] + adj
        }
    }
    
    # Correction Pass 2: Adjust for rounding errors
    diff_vars <- n_causal_sites - sum(chrom_summary$n_causal)
    if (diff_vars != 0) {
      if (diff_vars > 0) {
          eligible_chrs <- chrom_summary$n_causal < chrom_summary$available_sites
          pool <- chrom_summary$chr[eligible_chrs]
          pool_weights <- chrom_summary$total_weight[eligible_chrs]
      } else { # diff_vars < 0
          eligible_chrs <- chrom_summary$n_causal > 0
          pool <- chrom_summary$chr[eligible_chrs]
          pool_weights <- chrom_summary$total_weight[eligible_chrs]
      }
      if(length(pool) > 0) {
          adj_chr <- sample(pool, size=abs(diff_vars), replace=TRUE, prob=pool_weights)
          adj_table <- table(adj_chr)
          for (c in names(adj_table)) {
              chrom_summary$n_causal[chrom_summary$chr == c] <- chrom_summary$n_causal[chrom_summary$chr == c] + sign(diff_vars)*adj_table[c]
          }
      }
    }
    ##=====================================================================================================
    ### SECTION 4: SAMPLE CAUSAL SITES
    ##=====================================================================================================
    
    sampled_variants <- all_variants_annot %>%
           dplyr::left_join(dplyr::select(chrom_summary, chr, n_causal), by = "chr") %>%
           dplyr::group_by(chr) %>%
           dplyr::group_modify(~ {
               # .x is the data frame for the current group.
               # The n_causal value is constant for all rows in .x, so we take the first (can take any)
               num_to_sample <- .x$n_causal[1]
                   # Perform the sampling on the group's data frame (.x)
               dplyr::slice_sample(
                   .x,
                   n = num_to_sample,
                   weight_by = final_weight,
                   replace = FALSE
               )
           }) %>%
           dplyr::ungroup()

    causal_genotypes_list <- list()
    gds_groups <- split(sampled_variants, sampled_variants$gds_path)
    for (gds_path in names(gds_groups)) {
        gds_subset <- gds_groups[[gds_path]]
        if(is.null(gds_subset) || nrow(gds_subset) == 0) next
        cat(">>> Extracting causal sites:\n")
        dosage_mat <- safely_read_gds(gds_path, function(gds){
            SeqArray::seqSetFilter(gds, variant.id = gds_subset$variant.id.internal)
            # altDosage is samples x variants, so transpose to variants x samples
            t(SeqVarTools::altDosage(gdsobj = gds, use.names = TRUE, parallel = n_threads))
        })
        geno_df <- as.data.frame(dosage_mat)
        geno_df$chr <- as.character(gds_subset$chr)
        geno_df$pos <- as.integer(gds_subset$pos) 
        geno_df$ref <- as.character(gds_subset$ref)
        geno_df$alt <- as.character(gds_subset$alt)
        merged_df <- dplyr::inner_join(gds_subset, geno_df, by = c("chr", "pos", "ref", "alt"))
        causal_genotypes_list[[gds_path]] <- merged_df
    }

    ##=====================================================================================================
    ### SECTION 5: FINALIZE OUTPUT
    ##=====================================================================================================   
    ## Debugging lines:
    # cat("Number of rows in causal_genotypes_list:\n")
    # print(lapply(causal_genotypes_list, nrow))
    # cat("Number of columns in causal_genotypes_list:\n")
    # print(lapply(causal_genotypes_list, ncol))
   
    # causal_genotypes <- data.table::rbindlist(causal_genotypes_list, use.names = TRUE, fill = TRUE)
    causal_genotypes <- do.call("rbind", causal_genotypes_list)

    # Robustly derive sample and annotation columns from the final data frame
    all_cols <- names(causal_genotypes)
    annot_cols <- intersect(names(all_variants_annot), all_cols)
    sample_cols <- setdiff(all_cols, annot_cols)
    causal_genotypes <- causal_genotypes[, c(annot_cols, sample_cols)]
    
    causal_annot <- causal_genotypes[, annot_cols]

    output <- list(causal_annotation = causal_annot[,1:6], 
                    causal_genotypes = causal_genotypes[,-c(7:10)])
    return(output)
}

