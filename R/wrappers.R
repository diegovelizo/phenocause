#' Sample causal genetic variants (uniform mode)
#'
#' @description
#' Samples a specified number of causal variants from a set of GDS files. This function is a wrapper for the `uniform` mode of `sample_causal_sites()`, where variants are sampled uniformly at random. It can filter sites by missingness or minor allele count (MAC) before sampling and returns the genotypes and annotations of the selected variants.
#'
#' @details
#' This function triggers the `uniform` sampling mode, which is the default behavior. Variants are sampled uniformly at random from the pool of available sites after any filtering is applied.
#'
#' @param gds_paths A character vector of paths to the GDS files. When `length > 1`, each file is expected to contain data for a single chromosome.
#' @param n_causal_sites A single positive integer specifying the total number of causal sites to sample across all GDS files.
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
#' @seealso 
#' \code{\link{sample_causal_sites.ldak}}, \code{\link{sample_causal_sites.custom}}, \code{\link{sample_causal_sites.collider}}
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
#' ## Example: Uniform sampling of causal sites
#' ##================================================================
#' cat(">>> Running Example: Uniform Sampling\n")
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#' n_causal_total <- 200
#' uniform_sites <- sample_causal_sites.uniform(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_total,
#'   n_threads = 2
#' )
#' stopifnot(nrow(uniform_sites$causal_genotypes) == n_causal_total)
#' cat("OK: Uniform sampling returned the correct number of sites.\n")
#' }
sample_causal_sites.uniform <- function(gds_paths,
                                        n_causal_sites,
                                        missingness_threshold = NULL,
                                        mac_threshold = NULL,
                                        n_threads = 1) {
  
  # This is a simple wrapper that calls the main function, ensuring all
  # other mode-specific parameters are NULL.
  sample_causal_sites(
    gds_paths = gds_paths,
    n_causal_sites = n_causal_sites,
    missingness_threshold = missingness_threshold,
    mac_threshold = mac_threshold,
    n_threads = n_threads,
    ldak_ld_weights_paths = NULL,
    sampling_weights_paths = NULL,
    ldak_chrs = NULL,
    pca_matrix = NULL
  )
}

#' Sample causal genetic variants (ldak mode)
#'
#' @description
#' Samples a specified number of causal variants from a set of GDS files. This function is a wrapper for the `ldak_weights` mode of `sample_causal_sites()`, where variants are sampled using weights from LDAK software. It can filter sites by missingness or minor allele count (MAC) before sampling and returns the genotypes and annotations of the selected variants.
#'
#' @details
#' This function triggers the `ldak_weights` sampling mode. It robustly handles cases where the number of GDS files and weight files do not match, under the condition that one of the inputs is a single, consolidated file. If a matched list of per-chromosome GDS and 2-column LDAK files are provided, the chromosome for each weight file is inferred automatically.
#'
#' @param gds_paths A character vector of paths to the GDS files. When `length > 1`, each file is expected to contain data for a single chromosome.
#' @param n_causal_sites A single positive integer specifying the total number of causal sites to sample across all GDS files.
#' @param ldak_ld_weights_paths A character vector of paths to LDAK weights files. The file must have either two 
#'  or three columns. 
#' - Two-columns files: Raw (no header) *.short output file by `ldak --calc-weights`. The columns are `rsid` and `weight`.  
#' - Three-columns files: File with a header, containing columns `chr` `rsid` and `weight` (case-insensitive). 
#' See `ldak_chrs` for how chromosome information is handled for 2-column files. Cannot be used with `sampling_weights_paths`. 
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
#' @seealso 
#' \code{\link{sample_causal_sites.uniform}}, \code{\link{sample_causal_sites.custom}}, \code{\link{sample_causal_sites.collider}}
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
#' ## Example: Weighted sampling using LDAK weights
#' ##================================================================
#' cat(">>> Running Example: LDAK Weighted Sampling\n")
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#' data(phenocause.ldak_weights)
#' temp_dir <- tempdir()
#' 
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
#' # Since the number of GDS files and weight files match, the chromosome
#' # for each weight file can be inferred automatically.
#' ldak_sites <- sample_causal_sites.ldak(
#'   gds_paths = gds_paths,
#'   n_causal_sites = 300,
#'   ldak_ld_weights_paths = temp_ldak_paths,
#'   n_threads = 2
#' )
#' stopifnot(nrow(ldak_sites$causal_genotypes) == 300)
#' cat("OK: LDAK weighted sampling returned the correct number of sites.\n")
#' file.remove(temp_ldak_paths)
#' }
sample_causal_sites.ldak <- function(gds_paths,
                                     n_causal_sites,
                                     ldak_ld_weights_paths,
                                     ldak_chrs = NULL,
                                     stabilize_sampling_weights = TRUE,
                                     weights_power = 1,
                                     missingness_threshold = NULL,
                                     mac_threshold = NULL,
                                     n_threads = 1) {

  # This is a simple wrapper that calls the main function, ensuring all
  # other mode-specific parameters are NULL.
  sample_causal_sites(
    gds_paths = gds_paths,
    n_causal_sites = n_causal_sites,
    ldak_ld_weights_paths = ldak_ld_weights_paths,
    ldak_chrs = ldak_chrs,
    stabilize_sampling_weights = stabilize_sampling_weights,
    weights_power = weights_power,
    missingness_threshold = missingness_threshold,
    mac_threshold = mac_threshold,
    n_threads = n_threads,
    sampling_weights_paths = NULL,
    pca_matrix = NULL
  )
}


#' Sample causal genetic variants (custom weights mode)
#'
#' @description
#' Samples a specified number of causal variants from a set of GDS files. This function is a wrapper for the `custom_weights` mode of `sample_causal_sites()`, where variants are sampled using user-provided, custom weights. It can filter sites by missingness or minor allele count (MAC) before sampling and returns the genotypes and annotations of the selected variants.
#'
#' @details
#' This function triggers the `custom_weights` sampling mode. It robustly handles cases where the number of GDS files and weight files do not match, under the condition that one of the inputs is a single, consolidated file (e.g., one GDS file for all chromosomes, or one weight file for all chromosomes).
#'
#' The join between weights/annotations and genotypes is performed using a robust composite key (`chr`, `pos`, `ref`, `alt`) to ensure correctness.
#'
#' @param gds_paths A character vector of paths to the GDS files. When `length > 1`, each file is expected to contain data for a single chromosome.
#' @param n_causal_sites A single positive integer specifying the total number of causal sites to sample across all GDS files.
#' @param sampling_weights_paths A character vector of paths to custom sampling weight files. Files must have a header and contain mandatory columns (case-insensitive): `chr`, `pos`, and `weight`. The columns `rsid`, `ref`, and `alt` are optional; a warning is issued if they are missing.
#' @param stabilize_sampling_weights A logical flag. If `TRUE` (default), a small constant (10^-6) is added to all weights to prevent zero-probability sampling.
#' @param weights_power A numeric value used to raise the sampling weights to a power. Default is 1.
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
#' @seealso 
#' \code{\link{sample_causal_sites.uniform}}, \code{\link{sample_causal_sites.ldak}}, \code{\link{sample_causal_sites.collider}}
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
#' ## Example: Weighted sampling with a custom weights file
#' ##================================================================
#' cat(">>> Running Example: Custom Weights File\n")
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
#' custom_sites <- sample_causal_sites.custom(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_custom,
#'   sampling_weights_paths = tmp_custom_files
#' )
#'
#' stopifnot(nrow(custom_sites$causal_genotypes) == n_causal_custom)
#' cat("OK: Custom weighted sampling returned the correct number of sites.\n")
#' file.remove(tmp_custom_files)
#' }
sample_causal_sites.custom <- function(gds_paths,
                                       n_causal_sites,
                                       sampling_weights_paths,
                                       stabilize_sampling_weights = TRUE,
                                       weights_power = 1,
                                       missingness_threshold = NULL,
                                       mac_threshold = NULL,
                                       n_threads = 1) {
  
  # This is a simple wrapper that calls the main function, ensuring all
  # other mode-specific parameters are NULL.
  sample_causal_sites(
    gds_paths = gds_paths,
    n_causal_sites = n_causal_sites,
    sampling_weights_paths = sampling_weights_paths,
    stabilize_sampling_weights = stabilize_sampling_weights,
    weights_power = weights_power,
    missingness_threshold = missingness_threshold,
    mac_threshold = mac_threshold,
    n_threads = n_threads,
    ldak_ld_weights_paths = NULL,
    ldak_chrs = NULL,
    pca_matrix = NULL
  )
}


#' Sample causal genetic variants (collider mode)
#'
#' @description
#' Samples causal variants from a set of GDS files based on their correlation with principal components (PCs). This function is a wrapper for the `collider` mode of `sample_causal_sites()`. It is designed to generate sets of variants that can be used to study and benchmark methods for handling collider bias in genome-wide association studies.
#'
#' @details
#' This function triggers the `collider` sampling mode. It first calculates the squared multiple correlation (R^2) between each variant and the provided principal components. It then stratifies variants by minor allele frequency (MAF) and samples two parallel sets of causal variants:
#' \itemize{
#'   \item A **"structured" set** from variants with high R^2 values. Phenotypes generated using this set are expected to be susceptible to collider bias when PCs are included as covariates in a GWAS model.
#'   \item A **"non-structured" set** from variants with low R^2 values. This set serves as a negative control, as phenotypes generated from it should not be subject to the same bias.
#' }
#'
#' @param gds_paths A character vector of paths to the GDS files. When `length > 1`, each file is expected to contain data for a single chromosome.
#' @param n_causal_sites A single positive integer specifying the total number of causal sites to sample for *each* set (structured and non-structured).
#' @param pca_matrix A data frame or matrix containing principal components for all samples. The first column must be named `sample.id` (case-insensitive), and all subsequent columns must contain the PC scores (e.g., `PC1`, `PC2`, ...). Providing this argument activates the `collider` sampling mode.
#' @param collider_r2_candidate_threshold A numeric quantile between 0 and 1. Variants with an R^2 value above this quantile within their MAF stratum are placed in the "structured" set. Default is `0.975`.
#' @param collider_r2_control_threshold A numeric quantile between 0 and 1. Variants with an R^2 value below this quantile are placed in the "non-structured" set. Default is `0.025`.
#' @param missingness_threshold A numeric value between 0 and 1. Variants with a missingness rate greater than this threshold will be excluded before sampling. Default is `NULL`.
#' @param mac_threshold A positive integer. Variants with a minor allele count (MAC) less than this threshold will be excluded before sampling. Default is `NULL`.
#' @param n_threads A positive integer specifying the number of threads for parallel operations on GDS files. Default is 1.
#'
#' @return
#' A nested list containing the annotations and genotypes for the structured and non-structured variant sets.
#' \describe{
#'   \item{`annotations`}{A list of four data frames:
#'     \itemize{
#'       \item `causal_structured`: Annotations for the `n_causal_sites` variants sampled from the high-R^2 set.
#'       \item `causal_non_structured`: Annotations for the `n_causal_sites` variants sampled from the low-R^2 set.
#'       \item `non_causal_structured`: Annotations for the remaining high-R^2 variants not chosen as causal.
#'       \item `non_causal_non_structured`: Annotations for the remaining low-R^2 variants not chosen as causal.
#'     }
#'   }
#'   \item{`genotypes`}{A list of two data frames, providing the genotype matrices for both causal sets:
#'     \itemize{
#'       \item `causal_structured`: Genotypes for the `n_causal_sites` variants from the high-R^2 set.
#'       \item `causal_non_structured`: Genotypes for the `n_causal_sites` variants from the low-R^2 set.
#'     }
#'   }
#' }
#'
#' @references
#' Grinde, K.E., Browning, B.L., Reiner, A.P., Thornton, T.A. & Browning, S.R. (2024) Adjusting for principal components can induce collider bias in genome-wide association studies. *PLOS Genetics*, 20(12), e1011242.
#'
#' @seealso 
#' \code{\link{sample_causal_sites.uniform}}, \code{\link{sample_causal_sites.ldak}}, \code{\link{sample_causal_sites.custom}}
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
#' ## Example: Collider mode sampling
#' ##================================================================
#' cat(">>> Running Example: Collider Mode Sampling\n")
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#' 
#' # Goal: Select two parallel sets of causal variants: one "structured" set
#' # that is highly correlated with PCs, and one "non-structured" set that is not.
#'
#' # 1. Load metadata containing pre-computed PCs
#' data(phenocause.metadata)
#' pca_matrix <- phenocause.metadata[, c("sample.id", paste0("PC", 1:10))]
#' n_causal <- 50
#' 
#' # 2. Run sampling in collider mode
#' collider_sites <- sample_causal_sites.collider(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal,
#'   pca_matrix = pca_matrix,
#'   n_threads = 2
#' )
#'
#' # 3. Validate the selection by comparing R^2 summaries between the two sets
#' cat("\nSummary of R^2 for structured causal sites (high-R^2 set):\n")
#' summary(collider_sites$annotations$causal_structured$r2)
#'
#' cat("\nSummary of R^2 for non-structured causal sites (low-R^2 set):\n")
#' summary(collider_sites$annotations$causal_non_structured$r2)
#' 
#' # 4. Check that differences in R^2 are not due to differences in MAF
#' cat("\nSummary of MAF for structured causal sites:\n")
#' summary(collider_sites$annotations$causal_structured$maf)
#' 
#' cat("\nSummary of MAF for non-structured causal sites:\n")
#' summary(collider_sites$annotations$causal_non_structured$maf)
#' 
#' stopifnot(nrow(collider_sites$genotypes$causal_structured) == n_causal)
#' stopifnot(nrow(collider_sites$genotypes$causal_non_structured) == n_causal)
#' cat("\nOK: Collider mode returned the correct number of sites for both sets.\n")
#' }
sample_causal_sites.collider <- function(gds_paths,
                                         n_causal_sites,
                                         pca_matrix,
                                         collider_r2_candidate_threshold = 0.975,
                                         collider_r2_control_threshold = 0.025,
                                         missingness_threshold = NULL,
                                         mac_threshold = NULL,
                                         n_threads = 1) {

  # This is a simple wrapper that calls the main function, ensuring all
  # other mode-specific parameters are NULL.
  sample_causal_sites(
    gds_paths = gds_paths,
    n_causal_sites = n_causal_sites,
    pca_matrix = pca_matrix,
    collider_r2_candidate_threshold = collider_r2_candidate_threshold,
    collider_r2_control_threshold = collider_r2_control_threshold,
    missingness_threshold = missingness_threshold,
    mac_threshold = mac_threshold,
    n_threads = n_threads,
    ldak_ld_weights_paths = NULL,
    sampling_weights_paths = NULL,
    ldak_chrs = NULL
  )
}

