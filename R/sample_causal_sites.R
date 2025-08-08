#' Sample causal genetic variants
#'
#' @description
#' Samples a specified number of causal variants from a set of GDS files. The function operates in one of several modes: uniform sampling, weighted sampling based on LDAK or custom user-provided weights, or a special mode to select variants based on their correlation with principal components to study collider bias. It can filter sites by missingness or minor allele count (MAC) before sampling and returns the genotypes and annotations of the selected variants.
#'
#' @details
#' The sampling mode is determined automatically based on the provided arguments:
#' \itemize{
#'   \item \strong{uniform}: Default mode. Triggered if no weight files or PCA matrix are provided. Variants are sampled uniformly at random from the pool of available sites after filtering.
#'   \item \strong{ldak_weights}: Triggered by providing paths in `ldak_ld_weights_paths`. Variants are sampled using weights from LDAK software.
#'   \item \strong{custom_weights}: Triggered by providing paths in `sampling_weights_paths`. Variants are sampled according to user-defined weights.
#'   \item \strong{collider}: Triggered by providing a `pca_matrix`. This mode first calculates the squared multiple correlation (R^2) between each variant and the provided principal components (PCs). It then stratifies variants by minor allele frequency (MAF) and selects causal variants exclusively from either the pool of variants with high or low R^2 values within each MAF tranche. Including PCs in an linear model where the sites in the high R^2 set ("structured sites") are causal is expected to induce collider bias (See Grinde et al. ).  The low R^2 set ("non-structured sites") can serve as a benchmarking control.
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
#' @param pca_matrix A data frame or matrix containing principal components for all samples. The first column must be named `sample.id` (case-insensitive), and all subsequent columns must contain the PC scores (e.g., `PC1`, `PC2`, ...). Providing this argument activates the `collider` sampling mode.
#' @param collider_r2_candidate_threshold A numeric quantile between 0 and 1. In `collider` mode, variants with an R² value above this quantile within their MAF stratum are considered candidates for being causal. Default is `0.975`.
#' @param collider_r2_control_threshold A numeric quantile between 0 and 1. In `collider` mode, variants with an R² value below this quantile are selected as a negative control set. Default is `0.025`.
#'
#' @return
#' For `uniform`, `custom_weights`, and `ldak_weights` modes, a list containing two data frames:
#' \describe{
#'   \item{`causal_annotation`}{A data frame with annotations for each of the `n_causal_sites` sampled variants.}
#'   \item{`causal_genotypes`}{A data frame with `n_causal_sites` rows and `(k + N)` columns, where `k` is the number of annotation columns and `N` is the number of samples. It contains the genotypes and annotations for the causal variants.}
#' }
#' For `collider` mode, the return value is a list containing two elements, where the first element is a nested list:
#' \describe{
#'   \item{`annotations`}{A list containing three data frames of annotations: `causal_sites` (variants sampled to be causal), `non_selected_candidates` (other high-R² variants not chosen), and `control_sites` (low-R² variants).}
#'   \item{`causal_genotypes`}{The genotype data frame for only the `causal_sites`.}
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
#' ## Example 3: Weighted sampling using LDAK weights
#' ##================================================================
#' cat(">>> Running Example 3: LDAK Weights from R Object\n")
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
#' # Here, we infer chromosomes from the matched GDS file list
#' # so `ldak_chrs` is not needed.
#' ldak_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = 300,
#'   ldak_ld_weights_paths = temp_ldak_paths,
#'   n_threads = 2
#' )
#' stopifnot(nrow(ldak_sites$causal_genotypes) == 300)
#' cat("OK: LDAK weighted sampling returned the correct number of sites.\n")
#' file.remove(temp_ldak_paths)
#'
#' ##================================================================
#' ## Example 4: Collider mode sampling
#' ##================================================================
#' cat(">>> Running Example 4: Collider Mode Sampling\n")
#' # Goal: Select causal variants that are either highly or weakly correlated with PCs
#' # The set with high R^2 ("structured sites") are expected to induce collider bias.
#' # The set with low R^2 ("non structured sites") can serve as a benchmarking control.
#'
#' # 1. Load metadata containing pre-computed PCs
#' data(phenocause.metadata)
#' pca_matrix <- phenocause.metadata[, c("sample.id", paste0("PC", 1:10))]
#' n_causal <- 50
#' # 2. Run sampling in collider mode
#' collider_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal,
#'   mac_threshold=20,
#'   missingness_threshold = 0.02,
#'   pca_matrix = pca_matrix,
#'   collider_r2_control_threshold = 0.025, ## Sites with in the lowest 2.5% will be returned as controls (no collider bias)
#'   collider_r2_candidate_threshold = 0.975,  ## sites in the top 2.5% will be used as candidates for the causal sites.
#'   n_threads = 2
#' )
#'
#' # 3. Validate the selection by comparing R^2 summaries
#' cat("\nSummary of R^2 for structured causal sites  (high R^2 set):\n")
#' summary(collider_sites$annotations$causal_structured$r2)
#'
#' cat("\nSummary of R^2 for non-structured causal sites (low R^2 set):\n")
#' summary(collider_sites$annotations$causal_non_structured$r2)
#' 
#' # 4. Make sure differences in R^2 are not due to differences in MAF: 
#' cat("\nSummary of MAF for structured causal sites (high R^2 set):\n")
#' summary(collider_sites$annotations$causal_structured$maf)
#' 
#' cat("\nSummary of MAF for non-structured causal sites (low R^2 set):\n")
#' summary(collider_sites$annotations$causal_non_structured$maf)
#' 
#' stopifnot(nrow(collider_sites$genotypes$causal_structured) == n_causal || nrow(collider_sites$genotypes$causal_non_structured) == n_causal )
#' cat("\nOK: Collider mode sampling returned the correct number of sites and pools.\n")
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
                              n_threads = 1,
                              # New parameters for collider mode
                              pca_matrix = NULL,
                              collider_r2_candidate_threshold = 0.975,
                              collider_r2_control_threshold = 0.025
                              ) {

    # Helper function to safely open and close GDS files
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

    # Determine sampling mode
    sampling_mode <- "uniform"
    if (!is.null(sampling_weights_paths)) sampling_mode <- "custom_weights"
    if (!is.null(ldak_ld_weights_paths)) sampling_mode <- "ldak_weights"
    if (!is.null(pca_matrix)) sampling_mode <- "collider"
    
    # General file path validation
    weights_paths <- NULL
    if (sampling_mode %in% c("custom_weights", "ldak_weights")) {
        weights_paths <- if (sampling_mode == "custom_weights") sampling_weights_paths else ldak_ld_weights_paths
        if (sum(file.exists(weights_paths)) != length(weights_paths)) stop("!!! One or more weight files are missing.\n")
        if (any(duplicated(weights_paths))) stop("!!! Found duplicate paths in the weights vector.")
        
        n_gds <- length(gds_paths)
        n_weights <- length(weights_paths)
        if (n_gds != n_weights && n_gds != 1 && n_weights != 1) {
            stop(glue::glue("!!! Mismatch in file counts. GDS files ({n_gds}) and weight files ({n_weights}) do not match, and neither is of length 1."))
        }
    }

    # GDS content validation
    gds_chrom_list <- lapply(gds_paths, safely_read_gds, fun = function(gds) unique(SeqArray::seqGetData(gds, "chromosome")))
    if (length(gds_paths) > 1 && any(sapply(gds_chrom_list, length) > 1)) {
        stop("!!! When multiple GDS files are provided, each file must contain data for only one chromosome.")
    } else if (length(gds_paths) == 1 && length(gds_chrom_list[[1]]) > 1) {
        warning("Single GDS file contains multiple chromosomes. Sampling will proceed from all of them.")
    }

    gds_sample_ids <- safely_read_gds(gds_paths[1], fun = function(gds) SeqArray::seqGetData(gds, "sample.id"))
    if (length(gds_paths) > 1) {
        for (i in 2:length(gds_paths)) {
            current_ids <- safely_read_gds(gds_paths[i], fun = function(gds) SeqArray::seqGetData(gds, "sample.id"))
            if (!identical(gds_sample_ids, current_ids)) stop("!!! Sample IDs are not consistent across all provided GDS files.")
        }
    }

    ### Collide-specific validations 
    if (sampling_mode == "collider") {
        if (!is.data.frame(pca_matrix) && !is.matrix(pca_matrix)) stop("!!! `pca_matrix` must be a data frame or matrix.")
        if (ncol(pca_matrix) < 2) stop("!!! `pca_matrix` must have at least two columns (sample.id and one PC).")
        if (tolower(names(pca_matrix)[1]) != "sample.id") stop("!!! The first column of `pca_matrix` must be 'sample.id'.")
        if (!all(startsWith(tolower(names(pca_matrix)[-1]), "pc"))) stop("!!! All columns in `pca_matrix` after the first must start with 'PC'.")
        if (!is.numeric(collider_r2_candidate_threshold) || length(collider_r2_candidate_threshold) != 1 || collider_r2_candidate_threshold < 0 || collider_r2_candidate_threshold > 1) stop("!!! `collider_r2_candidate_threshold` must be a single numeric value between 0 and 1.")
        if (!is.numeric(collider_r2_control_threshold) || length(collider_r2_control_threshold) != 1 || collider_r2_control_threshold < 0 || collider_r2_control_threshold > 1) stop("!!! `collider_r2_control_threshold` must be a single numeric value between 0 and 1.")
        if (collider_r2_candidate_threshold <= collider_r2_control_threshold) stop("!!! `collider_r2_candidate_threshold` must be greater than `collider_r2_control_threshold`.")
        if (!is.null(ldak_ld_weights_paths) || !is.null(sampling_weights_paths)) stop("!!! 'collider' mode is mutually exclusive with 'ldak_weights' and 'custom_weights' modes.")
        if (weights_power != 1 || !stabilize_sampling_weights) warning("`weights_power` and `stabilize_sampling_weights` are ignored in 'collider' mode.")

        pca_sample_ids <- pca_matrix[[1]]
        if (!all(sort(gds_sample_ids) == sort(pca_sample_ids))) stop("!!! Sample IDs in GDS files and `pca_matrix` do not match.")
    }

    ##=====================================================================================================
    ### MAIN LOGIC BRANCHING
    ##=====================================================================================================

    if (sampling_mode %in% c("uniform", "custom_weights", "ldak_weights")) {
        ## START OF NON-COLLIDER BRANCH
        # SECTION 2: AGGREGATE AND VALIDATE WEIGHTS
        weights_df <- NULL
        if (sampling_mode != "uniform") {
            if (sampling_mode == "custom_weights") {
                res <- lapply(seq_along(weights_paths), function(i) {
                    p <- weights_paths[i]; df <- data.table::fread(p, header = TRUE); names(df) <- tolower(names(df))
                    if (!all(c("chr", "pos", "weight") %in% names(df))) stop(glue::glue("!!! Custom weight file {p} missing mandatory columns."))
                    list(data = df, missing = setdiff(c("rsid", "ref", "alt"), names(df)), path = p)
                })
                weights_list <- lapply(res, `[[`, "data"); missing_info <- Filter(function(r) length(r$missing) > 0, res)
                if (length(missing_info) > 0) {
                    msg <- paste(sapply(missing_info, function(m) glue::glue("File '{basename(m$path)}' missing: {paste(m$missing, collapse=', ')}")), collapse="\n")
                    warning(glue::glue("Optional columns were missing from custom weight files:\n{msg}"))
                }
                weights_df <- data.table::rbindlist(weights_list, use.names = TRUE, fill = TRUE)
            } else if (sampling_mode == "ldak_weights") {
                if (!is.null(ldak_chrs) && length(ldak_chrs) != length(ldak_ld_weights_paths)) stop("!!! `ldak_chrs` length must match `ldak_ld_weights_paths` length.")
                weights_list <- lapply(seq_along(ldak_ld_weights_paths), function(i) {
                    p <- ldak_ld_weights_paths[i]; header_df <- data.table::fread(p, nrows = 1, header=FALSE)
                    if (ncol(header_df) == 3 && all(sapply(header_df, is.character))) {
                        df <- data.table::fread(p, header = TRUE); names(df) <- tolower(names(df))
                        if (!all(c("chr", "rsid", "weight") %in% names(df))) stop(glue::glue("!!! 3-column LDAK file {p} has incorrect headers."))
                    } else {
                        if (length(ldak_ld_weights_paths) > 1 && length(gds_paths) > 1 && length(ldak_ld_weights_paths) == length(gds_paths)) {
                            df <- data.table::fread(p, header = FALSE); names(df) <- c("rsid", "weight"); df$chr <- gds_chrom_list[[i]]
                        } else {
                            if (is.null(ldak_chrs)) stop("!!! 2-column LDAK files provided, but chromosome info cannot be inferred. Use the `ldak_chrs` parameter.")
                            df <- data.table::fread(p, header = FALSE); names(df) <- c("rsid", "weight"); df$chr <- ldak_chrs[i]
                        }
                    }
                    return(df)
                })
                weights_df <- data.table::rbindlist(weights_list, use.names = TRUE, fill = TRUE)
            }
            initial_rows <- nrow(weights_df)
            essential_cols <- intersect(c("chr", "pos", "ref", "alt", "rsid", "weight"), names(weights_df))
            weights_df <- na.omit(weights_df, cols = essential_cols)
            if (initial_rows - nrow(weights_df) > 0) warning(glue::glue("{initial_rows - nrow(weights_df)} rows with NA removed from weights data."))
            if (!is.numeric(weights_df$weight)) stop("!!! 'weight' column must be numeric.")
            if (any(weights_df$weight <= 0)) stop("!!! All weights must be strictly positive (> 0).")
            weights_df$chr <- as.character(weights_df$chr)
            if ("pos" %in% names(weights_df)) weights_df$pos <- as.integer(weights_df$pos)
            if ("rsid" %in% names(weights_df)) weights_df$rsid <- as.character(weights_df$rsid)
        }

        # SECTION 3: DEFINE SAMPLING POOL AND SITES PER CHROMOSOME
        variant_list <- lapply(gds_paths, function(p) {
            safely_read_gds(p, function(gds) {
                n_var <- SeqArray::seqSummary(gds, verbose = FALSE)$num.variant; if (n_var == 0) return(NULL)
                keep <- rep(TRUE, n_var)
                if (!is.null(mac_threshold)) keep <- keep & (SeqArray::seqAlleleCount(gds, minor = TRUE) >= mac_threshold)
                if (!is.null(missingness_threshold)) keep <- keep & (SeqArray::seqMissing(gds, per.variant = TRUE) <= missingness_threshold)
                if (!any(keep)) return(NULL)
                SeqArray::seqSetFilter(gds, variant.sel = keep)
                data.frame(
                    chr = as.character(SeqArray::seqGetData(gds, "chromosome")), pos = as.integer(SeqArray::seqGetData(gds, "position")),
                    rsid = as.character(SeqArray::seqGetData(gds, "annotation/id")), ref = as.character(SeqArray::seqGetData(gds, "$ref")),
                    alt = as.character(SeqArray::seqGetData(gds, "$alt")), variant.id.internal = SeqArray::seqGetData(gds, "variant.id"),
                    gds_path = p, stringsAsFactors = FALSE
                )
            })
        })
        variant_list <- Filter(Negate(is.null), variant_list)
        if (length(variant_list) == 0) stop("!!! No variants available after MAC/missingness filters.")
        all_variants_annot <- data.table::rbindlist(variant_list, use.names = TRUE, fill = TRUE)

        if (sampling_mode != "uniform") {
            join_by_cols <- intersect(c("chr", "pos", "rsid", "ref", "alt"), names(weights_df))
            all_variants_annot <- dplyr::inner_join(all_variants_annot, weights_df, by = join_by_cols)
        } else {
            all_variants_annot$weight <- 1
        }
        
        if (nrow(all_variants_annot) < n_causal_sites) stop(glue::glue("!!! Cannot sample {n_causal_sites} sites. Only {nrow(all_variants_annot)} variants are available after filtering/merging."))
        
        all_variants_annot$final_weight <- (all_variants_annot$weight + ifelse(stabilize_sampling_weights, 1e-6, 0))^weights_power
        
        chrom_summary <- all_variants_annot %>% group_by(chr) %>% summarise(total_weight = sum(final_weight), available_sites = dplyr::n(), .groups = 'drop')
        chrom_summary$n_causal <- round(n_causal_sites * (chrom_summary$total_weight / sum(chrom_summary$total_weight)))
        
        excess <- sum(pmax(0, chrom_summary$n_causal - chrom_summary$available_sites))
        chrom_summary$n_causal <- pmin(chrom_summary$n_causal, chrom_summary$available_sites)
        if (excess > 0) {
            warning("! Number of available variants was lower than requested for at least one chromosome. Redistributing sites.")
            eligible_chrs <- chrom_summary$n_causal < chrom_summary$available_sites
            if(any(eligible_chrs)) {
                prop_weights <- chrom_summary$total_weight[eligible_chrs]
                adj <- round(excess * (prop_weights / sum(prop_weights)))
                chrom_summary$n_causal[eligible_chrs] <- chrom_summary$n_causal[eligible_chrs] + adj
            }
        }
        
        diff_vars <- n_causal_sites - sum(chrom_summary$n_causal)
        if (diff_vars != 0) {
            pool_info <- if (diff_vars > 0) {
                subset(chrom_summary, n_causal < available_sites)
            } else {
                subset(chrom_summary, n_causal > 0)
            }
            if(nrow(pool_info) > 0) {
                adj_chr <- sample(pool_info$chr, size=abs(diff_vars), replace=TRUE, prob=pool_info$total_weight)
                adj_table <- table(adj_chr)
                for (c in names(adj_table)) {
                    chrom_summary$n_causal[chrom_summary$chr == c] <- chrom_summary$n_causal[chrom_summary$chr == c] + sign(diff_vars)*adj_table[c]
                }
            }
        }
        
        # SECTION 4 and 5: SAMPLE SITES AND GENERATE OUTPUT
        sampled_variants <- all_variants_annot %>%
               dplyr::left_join(dplyr::select(chrom_summary, chr, n_causal), by = "chr") %>%
               dplyr::group_by(chr) %>%
               dplyr::group_modify(~ {
                   num_to_sample <- .x$n_causal[1]; if (is.na(num_to_sample)) num_to_sample <- 0
                   dplyr::slice_sample(.x, n = num_to_sample, weight_by = final_weight, replace = FALSE)
               }) %>% dplyr::ungroup()

        gds_groups <- split(sampled_variants, sampled_variants$gds_path)
        causal_genotypes_list <- lapply(names(gds_groups), function(gds_path) {
            gds_subset <- gds_groups[[gds_path]]; if(is.null(gds_subset) || nrow(gds_subset) == 0) return(NULL)
            dosage_mat <- safely_read_gds(gds_path, function(gds){
                SeqArray::seqSetFilter(gds, variant.id = gds_subset$variant.id.internal)
                t(SeqVarTools::altDosage(gdsobj = gds, use.names = TRUE, parallel = n_threads))
            })
            geno_df <- as.data.frame(dosage_mat)
            geno_df$chr <- as.character(gds_subset$chr); geno_df$pos <- as.integer(gds_subset$pos)
            geno_df$ref <- as.character(gds_subset$ref); geno_df$alt <- as.character(gds_subset$alt)
            dplyr::inner_join(gds_subset, geno_df, by = c("chr", "pos", "ref", "alt"))
        })
        
        causal_genotypes <- data.table::rbindlist(Filter(Negate(is.null), causal_genotypes_list), use.names = TRUE, fill = TRUE)
        all_cols <- names(causal_genotypes)
 
        annot_cols <- intersect(names(all_variants_annot), all_cols)
        sample_cols <- setdiff(all_cols, annot_cols)
        cols_to_extract <-  c(annot_cols, sample_cols)
        causal_genotypes <- as.data.frame(causal_genotypes)[, cols_to_extract] %>% 
                                dplyr::arrange(chr, pos)        
        causal_annot <- as.data.frame(causal_genotypes)[, annot_cols] %>% 
                            dplyr::arrange(chr, pos)
                
        return(list(causal_annotation = causal_annot[,1:6], causal_genotypes = causal_genotypes[,-c(7:10)]))
        ## END OF NON-COLLIDER BRANCH

    } else if (sampling_mode == "collider") {
        ## START OF COLLIDER BRANCH
        
        # Step 1: Data harmonization and R-squared pre-computation
        # Harmonize PCA matrix samples with GDS samples
        common_samples <- intersect(gds_sample_ids, pca_matrix[[1]])
        if(length(common_samples) < nrow(pca_matrix) || length(common_samples) < length(gds_sample_ids)){
            warning("Subsetting PCA matrix and GDS samples to a common set of individuals.")
        }
        pca_matrix_filt <- pca_matrix[pca_matrix[[1]] %in% common_samples, ]
        pc_scores <- as.matrix(pca_matrix_filt[, -1])
        rownames(pc_scores) <- pca_matrix_filt[[1]]
        pc_scores <- pc_scores[gds_sample_ids[gds_sample_ids %in% common_samples], ] # Match GDS order

        # Get minimal annotation for all passing variants first
        minimal_annot_list <- lapply(gds_paths, function(p) {
            safely_read_gds(p, function(gds) {
                SeqArray::seqSetFilter(gds, sample.id = common_samples) 
                keep <- rep(TRUE, SeqArray::seqSummary(gds, verbose=F)$num.variant)
                if (!is.null(mac_threshold)) keep <- keep & (SeqArray::seqAlleleCount(gds, minor = TRUE) >= mac_threshold)
                if (!is.null(missingness_threshold)) keep <- keep & (SeqArray::seqMissing(gds, per.variant = TRUE) <= missingness_threshold)
                if (!any(keep)) return(NULL)
                SeqArray::seqSetFilter(gds, variant.sel = keep)
                
                data.frame(
                    chr = as.character(SeqArray::seqGetData(gds, "chromosome")),
                    variant.id.internal = SeqArray::seqGetData(gds, "variant.id"),
                    maf = SeqArray::seqAlleleFreq(gds, minor = TRUE),
                    gds_path = p,
                    stringsAsFactors = FALSE
                )
            })
        })
        minimal_annot <- data.table::rbindlist(Filter(Negate(is.null), minimal_annot_list))
        
        if(nrow(minimal_annot) == 0) stop("No variants remained after filtering.")
        # Define the function for block processing
        calculate_r2_chunk <- function(geno_chunk) {
            cm <- colMeans(geno_chunk, na.rm=TRUE)
            for (j in seq_len(ncol(geno_chunk))) {
                geno_chunk[is.na(geno_chunk[, j]), j] <- cm[j]
            }
            geno_chunk_centered <- scale(geno_chunk, center = TRUE, scale = FALSE)
            cor_matrix <- cor(geno_chunk_centered, pc_scores, use = "pairwise.complete.obs")
            rowSums(cor_matrix^2, na.rm = TRUE)
        }

        # Run seqBlockApply on all GDS files to get R-squared values
        cat(glue(">>> Prefiltering finished. Regressing genotypes on principal components using {n_threads} threads..."), "\n")
        r2_list <- lapply(gds_paths, function(p) {
            safely_read_gds(p, function(gds) {
                # Ensure the same variant and sample filters are applied as for minimal_annot
                cat(glue('>>> Processing file with chromosome(s) {unique(seqGetData(gds, "chromosome"))}.'), "\n")
                gds_variants <- data.frame(minimal_annot)[minimal_annot$gds_path == p, "variant.id.internal"]
                cat(glue('>> Number of variants to be regressed: {length(gds_variants)}.'), "\n")
                SeqArray::seqSetFilter(gds, sample.id = common_samples, variant.id = gds_variants)
                # Pre-fetch variant IDs in the exact order they will be processed
                variant_ids_in_block_order <- SeqArray::seqGetData(gds, "variant.id")

                r2_vectors <- SeqArray::seqBlockApply(gds,
                                        var.name = "$dosage_alt", 
                                        FUN = calculate_r2_chunk, 
                                        margin="by.variant", 
                                        parallel=n_threads,
                                        as.is="list")
                chr_name <- as.character(SeqArray::seqGetData(gds, "chromosome"))
                data.frame(chr = chr_name, 
                           variant.id.internal = variant_ids_in_block_order, 
                           r2 = unlist(r2_vectors))
            })
        })
        cat(">>> Finished regressing genotypes on principal components.\n")

        r2_df <- data.table::rbindlist(r2_list)
        all_snp_data <- dplyr::inner_join(minimal_annot, r2_df, by = c("chr", "variant.id.internal"))
        # Step 2: MAF Stratification & Candidate/Control Selection
        maf_bins <- c(0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.5)
        all_snp_data$maf_tranche <- cut(all_snp_data$maf, breaks = maf_bins, include.lowest = TRUE, right=TRUE)
        tranches <- split(all_snp_data, all_snp_data$maf_tranche)

        # Handle sparse tranches by merging upwards
        min_proportion <- min(1 - collider_r2_candidate_threshold, 
                                collider_r2_control_threshold)
        min_sites_per_tranche <- ceiling(1 / min_proportion)

        tranche_names <- names(tranches)
        if (length(tranches) > 1) {
            for (i in (length(tranches) - 1):1) {
                if (nrow(tranches[[i]]) < min_sites_per_tranche) { 
                    # Use the original tranche name for the warning message
                    warning(glue("MAF tranche '{tranche_names[i]}' has only {nrow(tranches[[i]])} variants. Merging with next tranche."))
                    # Merge the sparse tranche (i) into the next one (i+1)
                    tranches[[i+1]] <- rbind(tranches[[i+1]], tranches[[i]])
                    # Mark the now-empty tranche for removal
                    tranches[[i]] <- NULL
                }
            }
        }
        
        # After the loop, remove all tranches that were marked for removal (set to NULL)
        tranches <- Filter(Negate(is.null), tranches)

        tranches <- Filter(Negate(is.null), tranches)
        selected_pools <- lapply(tranches, function(tranche_df) {
            r2_quantiles <- quantile(tranche_df$r2, 
                                        probs = c(collider_r2_control_threshold, collider_r2_candidate_threshold), 
                                        na.rm=TRUE)

            list(candidates = tranche_df[tranche_df$r2 >= r2_quantiles[2], ], 
                        controls = tranche_df[tranche_df$r2 <= r2_quantiles[1], ])
        })
        candidate_pool_df <- data.table::rbindlist(lapply(selected_pools, `[[`, "candidates")) %>% 
                                dplyr::arrange(chr, variant.id.internal)

        control_pool_df <- data.table::rbindlist(lapply(selected_pools, `[[`, "controls")) %>% 
                                dplyr::arrange(chr, variant.id.internal)
       
        if (nrow(candidate_pool_df) < n_causal_sites) {
            stop(glue("!!! Insufficient collider candidates found ({nrow(candidate_pool_df)}) to sample {n_causal_sites} causal sites."))
        }
        # Step 3: Final Causal Site Sampling
        set.seed(1)
        causal_structured_indices     <- sort(sample(1:nrow(candidate_pool_df), size = n_causal_sites, replace = FALSE))
        causal_non_structured_indices <- sort(sample(1:nrow(control_pool_df), size = n_causal_sites, replace = FALSE))

        causal_structured_sites_df <- candidate_pool_df[causal_structured_indices, ]
        causal_non_structured_sites_df <- control_pool_df[causal_non_structured_indices, ]

        non_causal_structured_df <- candidate_pool_df[-causal_structured_indices, ]
        non_causal_non_structured_df <- control_pool_df[-causal_non_structured_indices, ]

        # Step 4: Final Output Assembly
        all_selected_ids <- unique(c(causal_structured_sites_df$variant.id.internal, 
                                        causal_non_structured_sites_df$variant.id.internal, 
                                        non_causal_structured_df$variant.id.internal,
                                        non_causal_non_structured_df$variant.id.internal))
        
        full_annot_list <- lapply(gds_paths, function(p) {
            safely_read_gds(p, function(gds) {
                SeqArray::seqSetFilter(gds, variant.id = all_selected_ids)
                ## Note: Because we are filtering only by variant.id but not by chr at this step,
                # we might extract unwanted sites for now. However, we get rid of them later on
                # when merging back by chr and variant.d.internal 
                if (length(SeqArray::seqGetData(gds, "variant.id"))== 0) return(NULL)
                data.frame(
                    chr = as.character(SeqArray::seqGetData(gds, "chromosome")), ## added!!
                    pos = as.integer(SeqArray::seqGetData(gds, "position")), 
                    rsid = as.character(SeqArray::seqGetData(gds, "annotation/id")),
                    ref = as.character(SeqArray::seqGetData(gds, "$ref")), 
                    alt = as.character(SeqArray::seqGetData(gds, "$alt")),
                    variant.id.internal = SeqArray::seqGetData(gds, "variant.id"), 
                    stringsAsFactors = FALSE
                )
            })
        })

        full_annot <- data.table::rbindlist(Filter(Negate(is.null), full_annot_list))

        # Join full annotation back to the three sets
        final_cols <- c("chr", "pos", "rsid", "ref", "alt", "maf", "r2", "variant.id.internal", "gds_path") 

        causal_structured_annot <- dplyr::left_join(causal_structured_sites_df, 
                                                full_annot, by = c("chr", "variant.id.internal")) %>%
                                    dplyr::select(dplyr::all_of(final_cols))

        causal_non_structured_annot <- dplyr::left_join(causal_non_structured_sites_df, 
                                                full_annot, by = c("chr", "variant.id.internal")) %>%
                                    dplyr::select(dplyr::all_of(final_cols))

        non_causal_structured_annot <- dplyr::left_join(non_causal_structured_df, 
                                                full_annot, by = c("chr", "variant.id.internal")) %>%
                                    dplyr::select(dplyr::all_of(final_cols))

        non_causal_non_structured_annot <- dplyr::left_join(non_causal_non_structured_df, 
                                                full_annot, by = c("chr", "variant.id.internal")) %>%
                                    dplyr::select(dplyr::all_of(final_cols))


        ##############################################################
        # Causal genotypes matrix: structured sites
        gds_groups <- split(causal_structured_annot, causal_structured_annot$gds_path)
        
        structured_causal_genotypes_list <- lapply(names(gds_groups), function(gds_path) {
            gds_subset <- gds_groups[[gds_path]]
            if(is.null(gds_subset) || nrow(gds_subset) == 0) return(NULL)
            dosage_mat <- safely_read_gds(gds_path, function(gds){
                SeqArray::seqSetFilter(gds, variant.id = gds_subset$variant.id.internal, sample.id = common_samples)
                t(SeqVarTools::altDosage(gdsobj = gds, use.names = TRUE, parallel = n_threads))
            })
            geno_df <- as.data.frame(dosage_mat)
            geno_df$variant.id.internal <- gds_subset$variant.id.internal 
            dplyr::inner_join(gds_subset, geno_df, by = "variant.id.internal")
        })
 
        structured_causal_genotypes <- data.table::rbindlist(Filter(Negate(is.null), 
                                                structured_causal_genotypes_list), 
                                            use.names=TRUE, fill=TRUE)
            
        all_cols <- names(structured_causal_genotypes)
        annot_cols <- intersect(names(causal_structured_annot), all_cols)
        sample_cols <- setdiff(all_cols, annot_cols)
        structured_causal_genotypes <- as.data.frame(structured_causal_genotypes)[, c(annot_cols, sample_cols)]
        
        ##############################################################
        # Causal genotypes matrix: NON structured sites
        gds_groups <- split(causal_non_structured_annot, causal_non_structured_annot$gds_path)
        
        non_structured_causal_genotypes_list <- lapply(names(gds_groups), function(gds_path) {
            gds_subset <- gds_groups[[gds_path]]
            if(is.null(gds_subset) || nrow(gds_subset) == 0) return(NULL)
            dosage_mat <- safely_read_gds(gds_path, function(gds){
                SeqArray::seqSetFilter(gds, variant.id = gds_subset$variant.id.internal, sample.id = common_samples)
                t(SeqVarTools::altDosage(gdsobj = gds, use.names = TRUE, parallel = n_threads))
            })
            geno_df <- as.data.frame(dosage_mat)
            geno_df$variant.id.internal <- gds_subset$variant.id.internal 
            dplyr::inner_join(gds_subset, geno_df, by = "variant.id.internal")
        })
 
        non_structured_causal_genotypes <- data.table::rbindlist(Filter(Negate(is.null), 
                                                non_structured_causal_genotypes_list), 
                                            use.names=TRUE, fill=TRUE)
            
        all_cols <- names(non_structured_causal_genotypes)
        annot_cols <- intersect(names(causal_non_structured_annot), all_cols)
        sample_cols <- setdiff(all_cols, annot_cols)
        non_structured_causal_genotypes <- as.data.frame(non_structured_causal_genotypes)[, c(annot_cols, sample_cols)]
        ##############################################################

        annotations <- list( causal_structured = causal_structured_annot[,1:8], 
                            causal_non_structured = causal_non_structured_annot[,1:8], 
                            non_causal_structured = non_causal_structured_annot[,1:8], 
                            non_causal_non_structured = non_causal_non_structured_annot[,1:8])
        
        genotypes <- list( causal_structured = structured_causal_genotypes[,-9],
                           causal_non_structured = non_structured_causal_genotypes[,-9])

        return(list(annotations=annotations, genotypes=genotypes))
        ### END OF COLLIDER BRANCH
    }
}

