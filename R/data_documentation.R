#' Metadata for phenocause example genetic data
#'
#' @description
#' A data frame containing metadata for the 5,000 simulated individuals. This documentation also provides the complete generation details for the associated genetic data, which is stored externally as GDS files.
#'
#' @format
#' A `data.frame` with 5,000 rows and 22 columns:
#' \describe{
#'   \item{`sample.id`}{A unique identifier for each individual.}
#'   \item{`population`}{The simulated population label: Peruvian in Lima, Peru (PEL); Colombian in Medellin, Colombia (CLM); Mexican Ancestry in Los Angeles, California (MXL); Puerto Rican in Puerto Rico (PUR); Iberian in Spain (IBS); Yoruba in Ibadan, Nigeria (YRI); Han Chinese in Beijing, China (CHB); and individuals with Mexican ancestry from the 1000 Genomes Project (MXB).}
#'   \item{`PC1` to `PC20`}{The first 20 principal components of genetic variation.}
#' }
#'
#' @section Associated Genetic Data:
#' This metadata corresponds to the genetic data stored in three GDS files within the package. Their paths can be retrieved as follows:
#' `example_gds_paths <- system.file(paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause")`
#'
#' @details
#' The associated genetic data was generated using a detailed pipeline to simulate realistic patterns of admixed ancestry. The principal components in this metadata object were calculated from that genetic data.
#'
#' \strong{1. Simulation pipeline}
#'
#' The raw DNA sequences and genealogies were simulated using `msprime` (Kelleher et al., 2016; Baumdicker et al., 2022). The resulting tree sequence data was exported to VCF format. These VCFs were then processed using `bcftools` and `PLINK` (Chang et al., 2015) to create a final, analysis-ready dataset containing only biallelic SNPs that passed quality control (MAF >= 1% and LD pruning via `--indep-pairwise 500 50 0.2`). Finally, the pruned VCF for each chromosome was converted to the GDS format using `SeqArray::seqVCF2GDS()`.
#'
#' \strong{2. Msprime simulation parameters}
#'
#' \itemize{
#'   \item \strong{Sample composition}: A total of 5,000 individuals were simulated from the 8 populations listed in the `@format` section.
#'   \item \strong{Demographic model}: The simulation used the complex demographic model of admixed Latin American populations developed by Medina-Muñoz et al. (2023).
#'   \item \strong{Recombination maps}: The human recombination maps for chromosomes 20, 21, and 22 from genome build GRCh38 were used. To create a smaller example dataset, only the first one-third of each chromosome's genetic length (in cM) was simulated.
#'   \item \strong{Coalescent model}: The Discrete-Time Wright-Fisher (DTWF) model was used. This model was chosen over the standard Hudson model because it more accurately handles the large sample sizes and recent admixture events present in the demographic model, preventing known genealogical distortions (Bhaskar et al., 2014; Nelson et al., 2020).
#'   \item \strong{Mutation model}: A constant mutation rate of \eqn{1.25 \times 10^{-8}} per base pair per generation was used, with the Jukes-Cantor (JC69) nucleotide substitution model (Jukes & Cantor, 1969).
#' }
#'
#' @seealso
#' \code{\link[phenocause]{simulate_phenotype}}, \code{\link[phenocause]{sample_causal_sites}}
#'
#' @references
#' Baumdicker, F. et al. Efficient ancestry and mutation simulation with msprime 1.0. Genetics 220, iyab229 (2022).
#'
#' Bhaskar, A., Wang, Y. X. & Song, Y. S. Distortion of genealogical properties when the sample is very large. PNAS 111, 2361-2366 (2014).
#'
#' Chang, C. C. et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience 4, 7 (2015).
#'
#' Jukes, T. H. & Cantor, C. R. Evolution of Protein Molecules. in Mammalian protein metabolism vol. 3 21–132 (Academic Press, New York, 1969).
#'
#' Kelleher, J., Etheridge, A. M. & McVean, G. Efficient Coalescent Simulation and Genealogical Analysis for Large Sample Sizes. PLOS Comput. Biol. 12, e1004842 (2016).
#'
#' Medina-Muñoz, S. et al. Demographic modeling of admixed Latin American populations from whole genomes. Am J Hum Genet 110, 1503-1518 (2023).
#'
#' Nelson, D., Kelleher, J. & Ragsdale, A. P. Accounting for long-range correlations in genome-wide simulations of large cohorts. PLOS Genet 16, e1008619 (2020).
#'
#' @examples
#' # 1. Get paths to the associated GDS files
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"),
#'   package = "phenocause"
#' )
#'
#' # 2. Open one GDS file to extract sample IDs
#' gds <- SeqArray::seqOpen(gds_paths[1])
#' sample_ids_in_gds <- SeqArray::seqGetData(gds, "sample.id")
#' SeqArray::seqClose(gds)
#'
#' # 3. Load the metadata object
#' data(phenocause.metadata)
#'
#' # 4. Print a small section of the metadata
#' phenocause.metadata[1:6, 1:8]
#'
#' # 5. Test that the sample IDs in the metadata and GDS file match
#' all_match <- all(sort(phenocause.metadata$sample.id) == sort(sample_ids_in_gds))
#' cat("Sample IDs in metadata and GDS file match:", all_match, "\n")
#'
"phenocause.metadata"


#' SNP LD weights from the phenocause example dataset
#'
#' @description
#' A list of data frames containing LDAK-computed Linkage Disequilibrium (LD) weights for all SNPs in the example genetic dataset.
#'
#' @format
#' A list of three `data.frame` objects with two columns, corresponding to the *weights.short files generated by ldak `--calc-weights-all`:
#' \describe{
#'   \item{`rsid`}{The SNP identifier.}
#'   \item{`ldak_weight`}{The LDAK-computed LD weight.}
#' }
#'
#' @details
#' An LDAK weight for a given SNP measures the amount of local LD between that SNP and its neighbors. Within `phenocause`, these weights can be used to define sampling probabilities for causal sites, allowing for simulations where causal variants are enriched in either high- or low-LD regions.
#'
#' @seealso
#' \code{\link[phenocause]{sample_causal_sites}}, \code{\link[phenocause]{metadata}}
#'
#' @examples
#' # Load the weights data
#' data(phenocause.ldak_weights)
#' # View chromosome names
#' names(phenocause.ldak_weights)
#' # View weights for chromosome 1
#' head(phenocause.ldak_weights[["chr1"]])
#'
"phenocause.ldak_weights"

