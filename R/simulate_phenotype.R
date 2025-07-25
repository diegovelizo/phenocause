#' Simulate Complex Phenotypes with Controlled Confounding
#'
#' @description
#' Implements a simulation framework to generate phenotypes under various genetic architectures. The primary purpose is to formally model and generate data where the conflation between population structure and environmental confounding—a pervasive issue in human genetics research—can be explicitly controlled and studied. For this purpose, the function implements several distinct, additive confounding models.
#'
#' @details
#' The simulation is structured around a flexible, additive model. All confounder effects are optional and are added to a base genetic model.
#'
#' \strong{1. Simulation framework}
#'
#' The function first simulates a quantitative trait (\eqn{Y_{q}}) as a sum of a genetic component (\eqn{G}), one or more optional confounding effects (\eqn{P_{k}}), and a random residual error (\eqn{E}):
#' \deqn{Y_{q} = G + \sum_{k} P_{k} + E}
#' For binary outcomes, this quantitative trait is treated as an underlying, unobserved **liability**. A threshold is applied to this liability based on the desired prevalence to generate a binary case/control phenotype.
#'
#' \strong{2. The base genetic model}
#'
#' The genetic component \eqn{G} is the aggregate effect of \eqn{M} causal variants. The effect size \eqn{\beta_{j}} for variant \eqn{j} is drawn from a normal distribution where the variance is scaled by the variant's allele frequency \eqn{p_{j}} and a user-defined parameter, \eqn{\alpha}:
#' \deqn{\beta_{j} \sim \mathcal{N}(0, [2p_{j}(1-p_{j})]^{\alpha})}
#' The total genetic effect for an individual \eqn{i} is the sum of their dosage-weighted effects: \eqn{G_{i} = \sum_{j=1}^{M} g_{ij}\beta_{j}}. The residual error \eqn{E} is drawn from \eqn{\mathcal{N}(0, \sigma^{2}_{E})}, where the variance \eqn{\sigma^{2}_{E}} is calculated to satisfy the target heritability (\eqn{h^2}).
#'
#' \strong{3. Confounding models}
#'
#' The framework can simulate three distinct types of confounding effects that are added to the base model. They can be used in any combination.
#'
#' \strong{3.1. Categorical confounder}
#'
#' This models an observable categorical variable that influences the phenotype. The implementation distinguishes between nominal and ordinal confounders based on the data type of the input vector:
#' \itemize{
#'   \item \strong{Nominal (input is `character`)}: Each category is assigned an independent, random effect. This is suitable for variables without inherent order, such as country of birth.
#'   \item \strong{Ordinal (input is `factor`)}: The effects are simulated to be monotonically increasing across the defined factor levels. This is suitable for variables with a natural order, such as educational attainment (`"High School" < "Bachelors" < "PhD"`) or socioeconomic status.
#' }
#'
#' \strong{3.2. Non-specific genetic confounder}
#'
#' This models a quantitative confounder \eqn{X_{gc}} that is correlated with the total genetic effect \eqn{G}. Its effect on the phenotype is \eqn{P_{gc} = \gamma_{gc} X_{gc}}. The simulation is governed by:
#' \itemize{
#'   \item \code{genetic_confounder_coeff} (\eqn{\gamma_{gc}}): The linear effect of the confounder on the phenotype.
#'   \item \code{genetic_confounder_rel_var} (\eqn{w_{gc} = \frac{Var(P_{gc})}{Var(G)}}): The variance of the genetic confounding effect relative to the genetic variance.
#'   \item \code{genetic_confounder_cor} (\eqn{\rho_{gc} = Cor(G, X_{gc})}): The correlation between the genetic effect and the confounder.
#' }
#' The underlying confounder \eqn{X_{gc}} is simulated as \eqn{X_{gc} = b_{gc} G + E_{gc}}, where the coefficient \eqn{b_{gc} = \rho_{gc} \sqrt{w_{gc}} / |\gamma_{gc}|} is derived from the parameters.
#'
#' \strong{3.3. Ancestry-specific confounder}
#'
#' This models a quantitative confounder \eqn{X_{ac}} that is correlated with a specific axis of genetic ancestry \eqn{Q}. Its effect is \eqn{P_{ac} = \gamma_{ac} X_{ac}}. The simulation is governed by:
#' \itemize{
#'   \item \code{ancestry_component} (\eqn{Q}): The vector of ancestry component values for each individual.
#'   \item \code{ancestry_confounder_coeff} (\eqn{\gamma_{ac}}): The linear effect of the confounder on the phenotype.
#'   \item \code{ancestry_confounder_rel_var} (\eqn{w_{ac} = \frac{Var(P_{ac})}{Var(G)}}): The variance of the ancestry confounding effect relative to the genetic variance.
#'   \item \code{ancestry_confounder_cor} (\eqn{\rho_{ac} = Cor(Q, X_{ac})}): The correlation between the ancestry component and the confounder.
#' }
#' The underlying confounder \eqn{X_{ac}} is simulated as \eqn{X_{ac} = b_{ac} Q + E_{ac}}, where \eqn{b_{ac} = \frac{\rho_{ac} \sqrt{w_{ac}} SD(G)}{|\gamma_{ac}| SD(Q)}}.
#'
#' \strong{4. Liability threshold model for binary traits}
#'
#' For binary traits, a case/control status is determined by applying a threshold to the quantitative liability score \eqn{Y_q}. Two methods are available:
#' \itemize{
#'   \item \strong{`gaussian`}: Assumes the liability \eqn{Y_q} follows a normal distribution and uses \code{qnorm()} to find the threshold corresponding to the target prevalence. This is appropriate when the liability is expected to be unimodal and symmetric.
#'   \item \strong{`empirical`}: Uses the empirical cumulative distribution function (\code{quantile()}) of the simulated liability scores to find the threshold. This method is robust to non-normality and is **strongly recommended** when using categorical confounders, as they can introduce multimodality into the liability distribution.
#' }
#'
#' @param causal_genotypes A data frame of causal genotypes. The first \code{n_annot_cols} columns must contain variant annotations, followed by \eqn{N} columns for each sample's alternate allele dosage. By default, \code{n_annot_cols} is 6, which corresponds to the output by \code{phenocause::sample_causal_sites()}. 
#' @param heritability The target narrow-sense heritability (\eqn{h^2}) of the trait. A numeric value between 0 and 1.
#' @param alpha A numeric parameter modeling the relationship between MAF and variant effect size, as defined in the LDAK model. Default is -1, which corresponds to the GCTA model where all variants are assumed to explain the same amount of genetic variance.
#' @param phenotype The type of phenotype to simulate. Either `"quantitative"` (default) or `"binary"`.
#' @param prevalence The phenotype prevalence for binary traits. Required if `phenotype = "binary"`.
#' @param liab_dist The liability distribution model for binary traits. One of `"auto"` (default), `"gaussian"`, or `"empirical"`. If `"auto"`, it uses `"empirical"` when a categorical confounder is present, otherwise `"gaussian"`.
#' @param n_annot_cols The number of initial annotation columns in `causal_genotypes` to separate from the genotype matrix. Default is 5.
#' @param categorical_confounder A `character` (for nominal) or `factor` (for ordinal) vector representing a categorical confounder. Must have length \eqn{N}.
#' @param categorical_confounder_variance The variance explained by the categorical confounder's effect.
#' @param genetic_confounder_coeff The linear effect coefficient (\eqn{\gamma_{gc}}). All three `genetic_confounder_*` parameters must be provided together, or all must be `NULL`.
#' @param genetic_confounder_rel_var The variance of the genetic confounding effect relative to the genetic variance (\eqn{w_{gc}}).
#' @param genetic_confounder_cor The correlation (\eqn{\rho_{gc}}) between the genetic effect \eqn{G} and the confounder \eqn{X_{gc}}.
#' @param ancestry_component A numeric vector representing an ancestry component (\eqn{Q}). All four `ancestry_confounder_*` parameters must be provided together, or all must be `NULL`.
#' @param ancestry_confounder_coeff The linear effect coefficient (\eqn{\gamma_{ac}}).
#' @param ancestry_confounder_rel_var The variance of the ancestry confounding effect relative to the genetic variance (\eqn{w_{ac}}).
#' @param ancestry_confounder_cor The correlation (\eqn{\rho_{ac}}) between the ancestry component \eqn{Q} and the confounder \eqn{X_{ac}}.
#'
#' @return A list containing three data frames providing a complete record of the simulation:
#' \describe{
#'   \item{\code{phenotypes}}{A data frame with one row per sample, containing the final phenotypes and all their constituent components.}
#'   \item{\code{confounders}}{A data frame containing the simulated values of the underlying confounding variables.}
#'   \item{\code{coefficients}}{A data frame detailing the coefficients used to generate the confounding effects.}
#' }
#' @references
#' Speed, D., Hemani, G., Johnson, M.R. et al. Improved heritability estimation from genome-wide SNPs. Am J Hum Genet 91, 1011–1021 (2012).
#'
#' @importFrom dplyr %>% select inner_join
#' @importFrom stats rnorm var sd qnorm quantile model.matrix
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ##=====================================================
#' ## Example 1: Basic quantitative trait simulation
#' ##=====================================================
#' ## 1. We will first sample the causal sites from a GDS (Genomic Data Structure) file
#' # using the phenocause::sample_causal_sites() function.
#' # Note: To create a GDS file from a VCF file, see SeqArray::seqVCF2GDS 
#' # Get path  to the example GDS files:
#' gds_paths <- system.file(
#'   paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
#' )
#' set.seed(123)
#' # Sample the causal sites. The output will be with two data frames.
#' n_causal_sites <- 200 ## Define number of causal sites
#' uniform_sites <- sample_causal_sites(
#'   gds_paths = gds_paths,
#'   n_causal_sites = n_causal_sites,
#'   n_threads = 2
#' )
#'  causal_geno_df <- uniform_sites$causal_genotypes 
#'
#' ## 2. Simulate phenotype
#' target_h2 <- 0.6
#' sim_basic <- simulate_phenotype(
#'   causal_genotypes = causal_geno_df,
#'   heritability = target_h2,
#'   n_annot_cols = 6 ## Number of metadata columns in causal_geno_df
#' )
#'
#' # 3. Check output
#' observed_h2 <- var(sim_basic$phenotypes$g) / (var(sim_basic$phenotypes$g) + var(sim_basic$phenotypes$residual))
#' cat(sprintf("Target h2: %.3f | Observed h2: %.3f\n", target_h2, observed_h2))
#'
#' ##===========================================================
#' ## Example 2: Confounder correlated with polygenic score (g)
#' ##===========================================================
#' set.seed(456)
#' ## Goal: To simulate a trait caused by both genetic and non-genetic effects, where the non-genetic
#' # variable `X_gc` is correlated with the polygenic score `g`.
#' # Specifically, we will assume the model:
#' # Y = mu + g + P_gc + e =  mu + (X_causal_sites * beta_causal_sites) + (X_gc * genetic_confounder_coeff) + e
#' # Where Y is the phenotype, g=(X_causal_sites * beta_causal_sites) is the vector of total genetic effects
#' # And X_gc is a confounding variable correlated with g, such that:
#' # Cor(g, X_gc) = genetic_confounder_cor,
#' # And: Var(X_gc*genetic_confounder_coeff)/Var(g) = genetic_confounder_rel_var is the variance caused 
#' # by the confounder relative to the genetic variance.
#' # where  `genetic_confounder_cor`, `genetic_confounder_coeff` and `genetic_confounder_rel_var` are user-defined parameters.
#'  
#' ## 1. Define the parameters of the confounding effects
#' target_g_cor <- 0.7 ## Cor(g, X_gc) 
#' target_confounder_rel_var <- 0.2 
#'
#' sim_genetic_conf <- simulate_phenotype(
#'   causal_genotypes = causal_geno_df, ## reuse genotype data from Example 1
#'   heritability = 0.5,
#'   n_annot_cols = 6,
#'   genetic_confounder_coeff = 0.5,
#'   genetic_confounder_cor = target_g_cor,
#'   genetic_confounder_rel_var = target_confounder_rel_var
#' )
#'
#' # 2. Check output
#' g <- sim_genetic_conf$phenotypes$g
#' p_genetic_confounder <- sim_genetic_conf$phenotypes$p_genetic_confounder
#' x_genetic_confounder <- sim_genetic_conf$confounders$X_genetic_confounder ## Simulated confounder variable
#'
#' observed_g_cor <- cor(g, x_genetic_confounder)
#' observed_g_rel_var <- var(p_genetic_confounder) / var(g)
#'
#' cat(sprintf(
#'   "Target Cor(G, X_gc): %.3f | Observed: %.3f\n", target_g_cor, observed_g_cor
#' ))
#' cat(sprintf(
#'   "Target Var(P_gc)/Var(G): %.3f | Observed: %.3f\n", target_confounder_rel_var, observed_g_rel_var
#' ))
#'
#' ##====================================================================
#' ## Example 3: Confounder correlated with a specific ancestry component
#' ##====================================================================
#' set.seed(789)
#' ## Goal: To simulate a trait influenced by a confounding variable that
#' # is correlated with a specific ancestry componen (e.g. PC1, EUR % ancestry, etc.).
#' # The model is simular as in example 2, except that now the confounder is correlated 
#' # with a specific ancestry component (also see `Details`):
#' # Y = mu + g + P_ac + e =  mu + (X_causal_sites * beta_causal_sites) + (X_ac * ancestry_confounder_coeff) + e
#' # Use PCs from metadata:
#' data(phenocause.metadata) ## Load metadata
#' pc1 <- phenocause.metadata$PC1
#'
#' # 2. Simulate phenotype with an ancestry-correlated confounder
#' target_anc_cor <- 0.8 ## Cor(g, X_ac)
#' target_anc_rel_var <- 0.15 ## Var(P_ac)/Var(g)
#'
#' sim_anc_conf <- simulate_phenotype(
#'   causal_genotypes = causal_geno_df, ## reuse data from Example 1.
#'   heritability = 0.5,
#'   n_annot_cols = 6,
#'   ancestry_component = pc1,
#'   ancestry_confounder_coeff = 0.6,
#'   ancestry_confounder_cor = target_anc_cor,
#'   ancestry_confounder_rel_var = target_anc_rel_var
#' )
#'
#' # 3. Check output
#' g_anc <- sim_anc_conf$phenotypes$g
#' p_ancestry_confounder <- sim_anc_conf$phenotypes$p_ancestry_confounder
#' x_ancestry_confounder <- sim_anc_conf$confounders$X_ancestry_confounder
#'
#' observed_anc_cor <- cor(pc1, x_ancestry_confounder)
#' observed_anc_rel_var <- var(p_ancestry_confounder) / var(g_anc)
#'
#' cat(sprintf(
#'   "Target Cor(Q, X_ac): %.3f | Observed: %.3f\n", target_anc_cor, observed_anc_cor
#' ))
#' cat(sprintf(
#'   "Target Var(P_ac)/Var(G): %.3f | Observed: %.3f\n", target_anc_rel_var, observed_anc_rel_var
#' ))
#'
#'
#' ##==========================================================
#' ## Example 4: Binary trait with three different forms of
#' ## confounding (30 replicates)
#' ##==========================================================
#' library(ggplot2)
#' library(patchwork)
#' set.seed(1024)
#'
#' # 1. Define parameters 
#' target_prev <- 0.15 ## Prevalence of the binary trait.
#' n_replicates <- 30        ## Number of simulation replicates
#'
#' # Define a categorical confounder:
#' categorical_confounder <- as.character(phenocause.metadata$population) 
#' # Note: If `categorical_confounder` is of class "character", it will be
#' # treated as a nominal-scale variable. If it is of class factor,
#' # it will be treated as an ordinal-scale variable, with monotonically-increasing
#' # values for the categories as they appear in `levels(categorical_confounder)`.
#' # This can be useful to simulate confounding by, e.g., socio-economic status, 
#' # educational attainment, etc.
#'
#' # 2. Run simulation across replicates
#' rep_results <- replicate(n_replicates, {
#'   sim_full <- simulate_phenotype(
#'     causal_genotypes =  causal_geno_df, ## reuse data from example 1
#'     heritability = 0.4,
#'     phenotype = "binary",
#'     prevalence = target_prev,
#'     liab_dist = "empirical", # Recommended when using categorical confounders
#'     n_annot_cols = 6,
#'     # Confounder 1: Categorical
#'     categorical_confounder = categorical_confounder,
#'     categorical_confounder_variance = 0.05,
#'     # Confounder 2: Ancestry-correlated:
#'     ancestry_component = pc1,    ## Defined in example 3
#'     ancestry_confounder_coeff = 0.5,
#'     ancestry_confounder_rel_var = 0.1,
#'     ancestry_confounder_cor = 0.6,
#'     # Confounder 3: PGS-correlated:
#'     genetic_confounder_coeff = -0.4,
#'     genetic_confounder_rel_var = 0.1,
#'     genetic_confounder_cor = -0.5
#'   )
#'   # Return a named vector of key metrics for this replicate
#'   c(
#'     obs_prev = mean(sim_full$phenotypes$y_binary),
#'     obs_h2 = var(sim_full$phenotypes$g) / (var(sim_full$phenotypes$g) + var(sim_full$phenotypes$residual)),
#'     obs_anc_cor = cor(pc1, sim_full$confounders$X_ancestry),
#'     obs_anc_rel_var = var(sim_full$phenotypes$p_ancestry) / var(sim_full$phenotypes$g)
#'   )
#' }, simplify = FALSE)
#'
#' results_df <- as.data.frame(do.call(rbind, rep_results))
#'
#' # 3. Summarize and check results
#' print(summary(results_df))
#'
#' # 4. Plot distributions of observed metrics vs. targets
#' p1 <- ggplot(results_df, aes(y = obs_prev)) +
#'   geom_boxplot() +
#'   geom_hline(yintercept = target_prev, color = "red", linetype = "dashed") +
#'   labs(title = "Prevalence", y = "Observed", x = "") +
#'   coord_flip()
#'
#' p2 <- ggplot(results_df, aes(y = obs_h2)) +
#'   geom_boxplot() +
#'   geom_hline(yintercept = 0.4, color = "red", linetype = "dashed") +
#'   labs(title = "Heritability", y = "Observed", x = "") +
#'   coord_flip()
#'
#' p3 <- ggplot(results_df, aes(y = obs_anc_cor)) +
#'   geom_boxplot() +
#'   geom_hline(yintercept = 0.6, color = "red", linetype = "dashed") +
#'   labs(title = "Cor(Q, X_ac)", y = "Observed", x = "") +
#'   coord_flip()
#'
#' p4 <- ggplot(results_df, aes(y = obs_anc_rel_var)) +
#'   geom_boxplot() +
#'   geom_hline(yintercept = 0.1, color = "red", linetype = "dashed") +
#'   labs(title = "Var(P_ac)/Var(G)", y = "Observed", x = "") +
#'   coord_flip()
#'
#' # Combine plots
#' (p1 | p2) / (p3 | p4) + plot_annotation(title = "Simulation metrics across 30 replicates")
#' }
#' 

simulate_phenotype <- function(causal_genotypes,
                               heritability,
                               alpha = -1,
                               phenotype = "quantitative",
                               prevalence = NULL,
                               liab_dist = "auto",
                               n_annot_cols = 6,
                               ## Categorical confounder:
                               categorical_confounder = NULL,
                               categorical_confounder_variance = 0.1,
                               ## Non-specific genetic confounder:
                               genetic_confounder_coeff = NULL,
                               genetic_confounder_rel_var = NULL,
                               genetic_confounder_cor = NULL,
                               ## Ancestry-specific confounder:
                               ancestry_component = NULL,
                               ancestry_confounder_coeff = NULL,
                               ancestry_confounder_rel_var = NULL,
                               ancestry_confounder_cor = NULL
                               ) {


    ##=============================================================
    ### Section 1: Input validation
    ##=============================================================
    phenotype <- tolower(phenotype)
    liab_dist <- tolower(liab_dist)

    if (!phenotype %in% c("quantitative", "binary")) stop("!!! `phenotype` must be 'quantitative' or 'binary'.")
    if (phenotype == "binary") {
        if (is.null(prevalence)) stop("!!! `prevalence` must be specified for a binary phenotype.")
        if (!is.numeric(prevalence) || prevalence <= 0 || prevalence >= 1) stop("!!! `prevalence` must be a numeric value between 0 and 1.")
        if (!liab_dist %in% c("auto", "gaussian", "empirical")) stop("!!! `liab_dist` must be 'auto', 'gaussian', or 'empirical'.")
    }

    nongen_check <- sum(is.null(genetic_confounder_coeff), is.null(genetic_confounder_rel_var), is.null(genetic_confounder_cor))
    if (nongen_check > 0 && nongen_check < 3) stop("!!! Must specify all three `genetic_confounder_*` parameters or none.")

    anc_check <- sum(is.null(ancestry_component), is.null(ancestry_confounder_coeff), is.null(ancestry_confounder_rel_var), is.null(ancestry_confounder_cor))
    if (anc_check > 0 && anc_check < 4) stop("!!! Must specify all four `ancestry_confounder_*` parameters or none.")

    ##=============================================================
    ### Section 2: Input preprocessing & base simulation
    ##=============================================================
    genotype_matrix <- as.matrix(causal_genotypes[, -(1:n_annot_cols)])
    n_sites <- nrow(genotype_matrix)
    n_samples <- ncol(genotype_matrix)
    sample_ids <- colnames(genotype_matrix)

    P_categorical <- P_genetic <- P_ancestry <- NULL

    allele_freqs <- rowMeans(genotype_matrix, na.rm = TRUE) / 2
    V <- 2 * allele_freqs * (1 - allele_freqs)
    beta_causal <- stats::rnorm(n = n_sites, mean = 0, sd = sqrt(V^alpha))

    genotype_matrix <- t(genotype_matrix)
    genetic_effects <- as.numeric(
                                sweep(genotype_matrix, 2, colMeans(genotype_matrix, na.rm=TRUE), FUN ="-") %*%  beta_causal
                                )

    var_g <- stats::var(genetic_effects)
    var_e <- var_g * (1 - heritability) / heritability
    residuals <- stats::rnorm(n_samples, mean = 0, sd = sqrt(var_e))

    ##=============================================================
    ### Section 3: Confounder simulation
    ##=============================================================
    df_confounders <- data.frame(sample.id = sample_ids)
    df_coefficients <- data.frame(term = character(), coefficient = numeric())

    if (!is.null(categorical_confounder)) {
        n_levels <- length(unique(categorical_confounder))
        raw_coeffs <- stats::rnorm(n_levels, mean = 0, sd = 1)

        # If the confounder is a vector of class factor, assume ordinal scale
        # with monotonically-increasing coefficients across levels:
        if (is.factor(categorical_confounder)) {
            names(raw_coeffs) <- levels(categorical_confounder)
            raw_coeffs[levels(categorical_confounder)] <- sort(raw_coeffs)
            raw_coeffs <- unname(raw_coeffs) 
        }

        conf_matrix <- stats::model.matrix(~ categorical_confounder - 1) ## Use "-1" to avoid using a category as the reference
        # Calculate the raw effects and the scaling factor needed to achieve the target variance:
        raw_effects <- as.numeric(conf_matrix %*% raw_coeffs)
        scaling_factor  <- sqrt(categorical_confounder_variance / stats::var(raw_effects))
        adj_coeffs <- raw_coeffs * scaling_factor
        P_categorical <- as.numeric(conf_matrix %*% adj_coeffs)
        # Update the output data frames:
        df_confounders$X_categorical_confounder <- categorical_confounder
        df_coefficients <- rbind(df_coefficients, data.frame(term = colnames(conf_matrix), coefficient = adj_coeffs))
    }

    if (!is.null(genetic_confounder_coeff)) {
        b_gc <- genetic_confounder_cor * sqrt(genetic_confounder_rel_var) / abs(genetic_confounder_coeff)
        var_e_gc <- var_g * ((genetic_confounder_rel_var / genetic_confounder_coeff^2) - b_gc^2)
        if (var_e_gc < 0) stop("!!! Invalid genetic confounder parameters led to negative variance.")
        X_genetic_confounder <- b_gc * genetic_effects + stats::rnorm(n_samples, mean = 0, sd = sqrt(var_e_gc))
        P_genetic <- X_genetic_confounder * genetic_confounder_coeff
        df_confounders$X_genetic_confounder <- X_genetic_confounder
        df_coefficients <- rbind(df_coefficients, data.frame(term = "X_genetic_confounder", coefficient = genetic_confounder_coeff))
    }
    
    if (!is.null(ancestry_component)) {
        b_ac <- (ancestry_confounder_cor * sqrt(ancestry_confounder_rel_var) * stats::sd(genetic_effects)) / (abs(ancestry_confounder_coeff) * stats::sd(ancestry_component))
        var_e_ac <- (ancestry_confounder_rel_var * var_g / ancestry_confounder_coeff^2) - (b_ac^2 * stats::var(ancestry_component))
        if (var_e_ac < 0) stop("!!! Invalid ancestry confounder parameters led to negative variance.")
        X_ancestry_confounder <- b_ac * ancestry_component + stats::rnorm(n_samples, mean = 0, sd = sqrt(var_e_ac))
        P_ancestry <- X_ancestry_confounder * ancestry_confounder_coeff
        df_confounders$X_ancestry_confounder <- X_ancestry_confounder
        df_coefficients <- rbind(df_coefficients, data.frame(term = "X_ancestry_confounder", coefficient = ancestry_confounder_coeff))
    }

    ##=============================================================
    ### Section 4: Assemble final phenotype
    ##=============================================================
    y_quantitative <- genetic_effects + residuals
    if (!is.null(P_categorical)) y_quantitative <- y_quantitative + P_categorical
    if (!is.null(P_genetic)) y_quantitative <- y_quantitative + P_genetic
    if (!is.null(P_ancestry)) y_quantitative <- y_quantitative + P_ancestry

    df_phenotypes <- data.frame(
        sample.id = sample_ids,
        g = genetic_effects,
        p_categorical_confounder = if (!is.null(P_categorical)) P_categorical else NA,
        p_genetic_confounder = if (!is.null(P_genetic)) P_genetic else NA,
        p_ancestry_confounder = if (!is.null(P_ancestry)) P_ancestry else NA,
        residual = residuals,
        y_quantitative = y_quantitative
    )

    if (phenotype == "binary") {
        chosen_liab_dist <- "gaussian"
        if (liab_dist == "empirical") {
            chosen_liab_dist <- "empirical"
        } else if (liab_dist == "auto" && !is.null(categorical_confounder)) {
            chosen_liab_dist <- "empirical"
            warning("A categorical confounder was supplied with liab_dist='auto'. Setting liability distribution to 'empirical' to handle likely deviations from normality.")
        }
        
        if (chosen_liab_dist == "gaussian") {
            liability_threshold <- stats::qnorm(1 - prevalence, mean = mean(y_quantitative), sd = stats::sd(y_quantitative))
        } else { # empirical
            liability_threshold <- stats::quantile(y_quantitative, probs = 1 - prevalence)
        }
        df_phenotypes$y_binary <- as.integer(y_quantitative >= liability_threshold)
    }

    return(list(
        phenotypes = df_phenotypes,
        confounders = df_confounders,
        coefficients = df_coefficients
    ))
}


