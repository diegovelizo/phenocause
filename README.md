
# phenocause: A package for simulating complex phenotypes

A central challenge in genetic epidemiology is to disentangle the
genetic basis of a complex trait from the effects of environmental and
social factors, especially when those factors are correlated with
genetic ancestry or the true polygenic score. `phenocause` is a
simulation toolkit designed to address this problem by allowing
researchers to simulate phenotypes where the contributions of genetics
and different types of confounding can be
explicitly and precisely controlled. The package implements a two-step
workflow: first, sampling a set of causal variants to define a trait’s
genetic architecture, and second, simulating a final phenotype based on
those variants under a variety of user-specified confounding models.
This package includes example genotype data and metadata to allow users
to follow the tutorial.

## The simulation framework: Theoretical models

`phenocause` generates phenotypes using an additive model. All
confounding effects are optional and are added to a base genetic model.
The quantitative liability for a trait ($Y$) is simulated as the sum of
a genetic component ($g$), one or more optional confounding effects
($P_k$), and a random residual error ($E$).

The general model is:

$$Y = \mu + g + P_{cc} + P_{gc} + P_{ac} + E$$

Where:

- **$Y$** is the final quantitative phenotype, which is treated as an
  unobserved liability for binary (case/control) traits.  
- **$\mu$** is the baseline mean of the phenotype.  
- **$g$** is the individual’s true total genetic effect, or polygenic
  score, determined by their causal alleles.  
- **$P_{cc}$** is the additive effect from a **c**ategorical
  **c**onfounder (e.g., socioeconomic status).  
- **$P_{gc}$** is the additive effect from a quantitative confounder
  correlated with the **g**enetic component $g$.  
- **$P_{ac}$** is the additive effect from a quantitative confounder
  correlated with a specific **a**ncestry **c**omponent $Q$.  
- **$E$** is a normally-distributed residual error term, the variance of
  which is scaled to meet the target heritability.

### The genetic component ($g$)

The total genetic effect for an individual $i$ is the sum of their
dosages at $M$ causal variants, weighted by each variant’s effect size,
$\beta_j$:

$$g_i = \sum_{j=1}^{M} g_{ij}\beta_{j}$$

The effect sizes are drawn from a normal distribution whose variance is
dependent on the allele frequency $p_j$ of the variant, controlled by
the parameter $\alpha$. This model is adopted from the LDAK software.

$$\beta_{j} \sim \mathcal{N}(0, [2p_{j}(1-p_{j})]^{\alpha})$$

The $\alpha$ parameter models the relationship between allele frequency
and the variance explained by a variant. A value of $\alpha = -1$, the
default, corresponds to the standard GCTA model where all variants are
assumed to contribute equally to heritability, regardless of their
frequency. Other values can be used to model architectures where, for
example, low-frequency variants have larger effects.

### The confounding components ($P$)

`phenocause` can simulate three distinct types of confounding.

#### Categorical confounder ($P_{cc}$)

This models an observable categorical variable, $X_{cc}$, that
influences the phenotype. The total effect is
$P_{cc} = X_{cc}\beta_{cc}$. The simulation distinguishes between two
scales:

- **Nominal Scale**: When the input is a `character` vector (e.g.,
  country of birth, population label), each category is assigned an
  independent, random effect.

- **Ordinal Scale**: When the input is an ordered `factor` (e.g.,
  educational attainment levels), the effects are simulated to be
  monotonically increasing across the factor levels.

The total variance explained by this component is controlled by the
`categorical_confounder_variance` parameter.

#### Quantitative confounder correlated with the polygenic score ($P_{gc}$)

This models a scenario where a latent quantitative confounder, $X_{gc}$,
is correlated with the true polygenic score ($g$). This model serves as
an alternative hypothesis to the scepecific ancestry-correlated
confounder (see below). The simulation is governed by three parameters:

- `genetic_confounder_coeff` ($\gamma_{gc}$): The direct linear effect
  of the confounder on the phenotype, where
  $P_{gc} = X_{gc}\gamma_{gc}$.  
- `genetic_confounder_cor` ($\rho_{gc}$): The target correlation between
  the genetic effect and the confounder, $Cor(g, X_{gc})$.  
- `genetic_confounder_rel_var` ($w_{gc}$): The target variance of the
  confounding effect *relative to* the genetic variance,
  $Var(P_{gc})/Var(g)$.

#### Quantitative confounder correlated with a specific ancestry component ($P_{ac}$)

This models the classic confounding case in genetic epidemiology where a
latent quantitative confounder, $X_{ac}$, is correlated with a specific
axis of genetic ancestry, $Q$. $Q$ can be a principal component, an
ancestry proportion, or any other continuous measure of ancestry. This
is designed to simulate how social, cultural, or environmental exposures
that are correlated with ancestry can bias genetic association studies.
The simulation is governed by an analogous set of four parameters:

- `ancestry_component` ($Q$): The vector of ancestry values for each
  individual.  
- `ancestry_confounder_coeff` ($\gamma_{ac}$): The direct effect of the
  confounder, where $P_{ac} = X_{ac}\gamma_{ac}$.  
- `ancestry_confounder_cor` ($\rho_{ac}$): The target correlation
  between the ancestry component and the confounder, $Cor(Q, X_{ac})$.  
- `ancestry_confounder_rel_var` ($w_{ac}$): The target variance of the
  confounding effect relative to the genetic variance,
  $Var(P_{ac})/Var(g)$.

### Simulating binary traits

For binary (case/control) outcomes, the quantitative trait $Y$ is
treated as an unobserved liability. A threshold is applied to this
liability distribution, and individuals with a liability exceeding the
threshold are assigned “case” status. The threshold value is determined
by the user-specified `prevalence` of the trait.

## Installation

You can install `phenocause` from GitHub using `devtools`.

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("diegovelizo/phenocause")
```

## The two-step simulation workflow

The `phenocause` simulation pipeline consists of two main functions:
`sample_causal_sites()` followed by `simulate_phenotype()`. The first
function requires genotype data to be stored in the GDS (Genomic Data
Structure) format. (See the section below on “[Working with the Genomic
Data Structure (GDS)
Format](#working-with-the-genomic-data-structure-gds-format)” for
details on why this format is used and how to create your own files).

[Skip explanation and proceed to the step-by-step
tutorial](#step-by-step-tutorial)

## Working with the Genomic Data Structure (GDS) format

### Why GDS?

The package is designed to operate on genotype data stored in GDS files.
GDS files store the same informat as a VCF file; however, the GDS format
is optimized for high-performance computation and efficient access. Its
primary advantages are:

- **Efficiency:** Data is stored in on-disk chunks, allowing for much
  faster access to specific variants or samples compared to parsing a
  text-based VCF file.

- **Memory Management:** The entire dataset does not need to be loaded
  into memory. Instead, a connection is opened to the file on disk, and
  only the requested data slices are read.

- **On-Disk Filtering:** Filters (e.g., by minor allele count,
  missingness or sample ID) can be applied to the file connection
  itself, before loading the data, ensuring that subsequent operations
  only process the subset of data that meets the specified criteria.

### Obtaining GDS files

This package includes example GDS files that can be used to follow the
tutorial. The paths on your system can be retrieved with the following
command:

``` r
gds_paths <- system.file(
  paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
)
# print(gds_paths)
```

To use your own data, you must first convert it from VCF format to GDS
using the `seqVCF2GDS()` function from the `SeqArray` package.

``` r
library(SeqArray)

# Define path to an existing VCF file and the desired output path for the GDS file
vcf_path <- "/path/to/your/data.vcf.gz" 
gds_path <- "/path/to/your/output.gds" 

# Convert the VCF to GDS. GT (genotype) is the minimum required format.
seqVCF2GDS(vcf_path, gds_path, fmt.import="GT", storage.option="LZMA_RA", parallel = 4)
```

### GDS access demonstration \[optional\]

Once you have a GDS file, you can open a connection and perform various
high-performance operations.  
Note: This is just a demonstration of how GDS files can be efficiently
accessed. These operations are done under the hood by `phenocause`, so
you don’t need to deal with them.

``` r
library(SeqArray)
library(SeqVarTools)

# Use the first example GDS file
gds_path <- gds_paths[1] 

# 1. Open a connection to the GDS file (does not load data into memory)
gds <- seqOpen(gds_path)

# 2. Get basic information
sample_ids <- seqGetData(gds, "sample.id")
variant_ids <- seqGetData(gds, "variant.id") # These are internal indices
cat("First 6 sample IDs:", head(sample_ids), "\n")
#> First 6 sample IDs: tsk.1 tsk.2 tsk.3 tsk.4 tsk.5 tsk.6

# 3. Get a detailed summary of the file contents
# seqSummary(gds, "genotype")

# 4. Apply filters without loading data
# Get per-variant statistics
mac <- seqAlleleCount(gds, minor = TRUE)
missingness <- seqMissing(gds, per.variant = TRUE)

# Create boolean vectors for filtering
mac_filter <- mac >= 20
missingness_filter <- missingness <= 0.02
final_filter <- mac_filter & missingness_filter

# Apply the filter to the GDS connection. 
# Subsequent operations will only use variants that pass this filter.
seqSetFilter(gds, variant.sel = final_filter)
#> # of selected variants: 6,559

# 5. Access data only from the filtered sites
# Get annotations for the filtered variants
sites_annotation <-  data.frame(
    chr = seqGetData(gds, "chromosome"),
    pos = seqGetData(gds, "position"),
    rsid = seqGetData(gds, "annotation/id")
)
cat("\nAnnotations for the first 6 variants passing filters:\n")
#> 
#> Annotations for the first 6 variants passing filters:
print(head(sites_annotation))
#>   chr  pos rsid
#> 1   1 1400   20
#> 2   1 2337   40
#> 3   1 3028   58
#> 4   1 4259   88
#> 5   1 4432   90
#> 6   1 6100  128

# Get genotype dosages for the filtered variants
# The altDosage function returns a matrix of samples x variants
genotype_data <- altDosage(gdsobj = gds, use.names = TRUE, parallel = 2)
cat("\nDimensions of genotype matrix from filtered sites (samples x variants):\n")
#> 
#> Dimensions of genotype matrix from filtered sites (samples x variants):
print(dim(genotype_data))
#> [1] 5000 6559


# 6. Close the GDS file connection
seqClose(gds)
```

## Data included with `phenocause`

The package includes three data objects to support the tutorial. They
can be loaded with `data()`.

### Genetic data

The example GDS files contain genotypes for 5,000 individuals across
portions of chromosomes 20, 21, and 22. The data was simulated using
`msprime` with a demographic model of Latin American and reference
continental populations resembling the 1000 Genomes Project populations,
to provide a realistic background of population structure and admixture.

### Metadata (`phenocause.metadata`)

This data frame contains sample-level information for the 5,000
simulated individuals, including population labels and the first 20
principal components calculated from the genetic data.

``` r
data(phenocause.metadata)
cat("Metadata for the first 6 samples:\n")
#> Metadata for the first 6 samples:
print(head(phenocause.metadata[, 1:8]))
#>    sample.id population         PC1         PC2         PC3          PC4
#>       <char>     <char>       <num>       <num>       <num>        <num>
#> 1:  sample.1        PEL -0.00773583 -0.00267714 -0.00841580 -1.44793e-02
#> 2:  sample.2        PEL -0.00863668 -0.00816183 -0.00813049 -6.31525e-03
#> 3:  sample.3        PEL -0.00862303 -0.00937326 -0.00993566  3.54648e-02
#> 4:  sample.4        PEL -0.00817220 -0.00699110 -0.00997924  3.06836e-02
#> 5:  sample.5        PEL -0.00960193 -0.01718230 -0.01265610  5.95974e-05
#> 6:  sample.6        PEL -0.00630836 -0.01232940 -0.01014900 -1.62568e-02
#>            PC5         PC6
#>          <num>       <num>
#> 1:  0.04927020  0.01048540
#> 2: -0.03133970 -0.00317630
#> 3: -0.00470034 -0.00986719
#> 4: -0.01424360 -0.00314133
#> 5: -0.00729436  0.00823931
#> 6: -0.02224340  0.00748102
```

### LDAK weights (`phenocause.ldak_weights`)

This is a list of data frames containing pre-computed Linkage
Disequilibrium (LD) weights from the LDAK software. Each data frame
corresponds to a chromosome in the example dataset. These weights
measure the extent of local LD surrounding each variant and can be used
to simulate genetic architectures where causal variant probability is
related to LD patterns.

``` r
data(phenocause.ldak_weights)
cat("LDAK weights for the first 6 variants on chromosome 1:\n")
#> LDAK weights for the first 6 variants on chromosome 1:
# The data object is a list named by chromosome
print(head(phenocause.ldak_weights[["chr1"]]))
#>     rsid   weight
#>    <int>    <num>
#> 1:    20 0.888681
#> 2:    40 0.702670
#> 3:    58 0.838932
#> 4:    88 0.809851
#> 5:    90 0.632627
#> 6:   128 0.587576
```

## Step-by-step tutorial

This tutorial demonstrates the two-step workflow for simulating
phenotypes with `phenocause`. First, we select causal variants to define
the genetic architecture. Second, we use those variants to simulate a
phenotype with specified genetic and confounding effects.

Before starting, we load the required libraries and get the paths to the
example GDS files included with the package.

``` r
library(phenocause)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glue)
library(SeqArray)
library(SeqVarTools)
library(data.table)
```

``` r
gds_paths <- system.file(
  paste0("extdata/example.chr", 1:3, ".gds"), package = "phenocause"
)
```

### Part 1: `sample_causal_sites()` - Defining the genetic architecture

**Rationale:** The first step is to define the degree of polygenicity
and how the causal sites are distributed along the genome. For instance,
we can sample causal sites with uniform probability along the genome, or
we can oversample causal sites at high- or low-LD regions (e.g. using LD
weight files calculated with `ldak --calc-weights`), or at any other
region by providing custom weight files. The `sample_causal_sites()`
function is designed for this purpose, providing several methods to
reflect different hypotheses about how causal variants are distributed
across the genome.

#### Example 1.1: Uniform sampling (basic model)

**Rationale:** The simplest model of genetic architecture is one where
every variant has an equal *a priori* probability of being causal. This
is achieved by sampling with uniform probability along the genome. We
will sample 1,000 sites and use this set for the phenotype simulation
examples in Part 2. Note:

``` r
set.seed(123)
n_causal_sites <- 1000 

uniform_sites <- sample_causal_sites(
  gds_paths = gds_paths,
  n_causal_sites = n_causal_sites,
  n_threads = 2
)
#> # of selected variants: 6,559
#> # of selected variants: 3,851
#> # of selected variants: 4,208
#> >>> Extracting causal sites:
#> # of selected variants: 449
#> >>> Extracting causal sites:
#> # of selected variants: 263
#> >>> Extracting causal sites:
#> # of selected variants: 288

# The output is a list containing two data frames. We will use 'causal_genotypes'
# as the primary input for the phenotype simulation function.
causal_geno_df <- uniform_sites$causal_genotypes
cat("Dimensions of the sampled causal genotype matrix:\n")
#> Dimensions of the sampled causal genotype matrix:
print(dim(causal_geno_df))
#> [1] 1000 5006
```

#### Example 1.2: Weighted sampling (LD weights)

**Rationale:** Complex trait genetics are shaped by evolutionary forces
like negative selection, which may cause causal variants to be
preferentially located in certain regions of the genome (e.g., regions
of lower linkage disequilibrium). The LDAK model provides a framework
that helps us simulate such architectures. By using weights derived from
local LD calculated via `ldak --calc-weights`, we can enrich our sample
for causal variants in high-LD or low-LD regions. The raw weights
(\*.weight.short) files output by LDAK measure the LD weights with
values ranging from 0 for variants in high LD with its neighbors to 1
for variants in low LD with its neighbors. The user can further
fine-tune the level of enrichment by using the `weights_power` parameter
in the phenocause::sample_causal_sites() function. For example, setting
`weights_power = -0.25` up-weights variants with smaller LD weights,
thus enriching the sample for variants in regions of higher LD.
Conversely, using `weights_power = 0.25` will oversample causal sites at
low-LD regions.

**Input format:** The function is flexible regarding the format of LDAK
weight files. \* **3-Column with header:** A file with `chr`, `rsid`,
and `weight` columns. This is the most robust format. \* **2-Column raw
output:** The default `.weights.short` file from LDAK containing `rsid`
and `weight` columns with no header. The package handles this format in
two ways: 1. If the number of weight files matches the number of GDS
files, chromosome information is automatically inferred from the
corresponding GDS file. 2. If the file counts do not match, the
`ldak_chrs` argument is **required** to manually provide the chromosome
for each weight file.

``` r
# Load the example LDAK weights data shipped with the package
data(phenocause.ldak_weights)
temp_dir <- tempdir()

# Create temporary 2-column LDAK weight files (no header) to demonstrate
temp_ldak_paths <- sapply(names(phenocause.ldak_weights), function(chr_name) {
  file_path <- file.path(temp_dir, paste0(chr_name, ".weights.short"))
  fwrite(
    phenocause.ldak_weights[[chr_name]],
    file = file_path, sep = " ", col.names = FALSE
  )
  return(file_path)
})

# Since the number of weight files (3) matches the number of GDS files (3),
# the `ldak_chrs` parameter is not strictly necessary here, but we include it
# for demonstration.
ldak_sites <- sample_causal_sites(
  gds_paths = gds_paths,
  n_causal_sites = 500,
  ldak_ld_weights_paths = temp_ldak_paths,
  ldak_chrs = c(1, 2, 3),
  weights_power = -0.25, # Example to enrich for high-LD regions
  n_threads = 2
)
#> # of selected variants: 6,559
#> # of selected variants: 3,851
#> # of selected variants: 4,208
#> >>> Extracting causal sites:
#> # of selected variants: 225
#> >>> Extracting causal sites:
#> # of selected variants: 132
#> >>> Extracting causal sites:
#> # of selected variants: 143

cat("Number of sites sampled with LDAK weights:", nrow(ldak_sites$causal_genotypes), "\n")
#> Number of sites sampled with LDAK weights: 500
```

#### Example 1.3: Weighted sampling (custom weights)

**Rationale:** This mode offers maximum flexibility, allowing a
researcher to test hypotheses about any user-defined genomic feature.
For example, one could use weights based on functional annotation scores
(e.g., CADD, PhyloP) to simulate a trait where causal variants are
enriched in functionally important regions of the genome. The input file
must contain `chr`, `pos`, and `weight` columns.

``` r
# Create temporary custom weight files with random weights for demonstration
tmp_custom_files <- sapply(1:3, function(chr_idx) {
  tmp_file <- tempfile(fileext = ".tsv")
  gds <- seqOpen(gds_paths[chr_idx])
  on.exit(seqClose(gds))
  
  annot_df <- data.frame(
      chr = seqGetData(gds, "chromosome"),
      pos = seqGetData(gds, "position"),
      rsid = seqGetData(gds, "annotation/id"),
      ref = seqGetData(gds, "$ref"),
      alt = seqGetData(gds, "$alt")
  )
  # For this example, we assign random weights
  annot_df$weight <- abs(rnorm(n = nrow(annot_df), mean = 0, sd = 3))
  fwrite(annot_df, file = tmp_file, sep = "\t")
  return(tmp_file)
})

custom_sites <- sample_causal_sites(
  gds_paths = gds_paths,
  n_causal_sites = 500,
  sampling_weights_paths = tmp_custom_files
)
#> # of selected variants: 6,559
#> # of selected variants: 3,851
#> # of selected variants: 4,208
#> >>> Extracting causal sites:
#> # of selected variants: 222
#> >>> Extracting causal sites:
#> # of selected variants: 133
#> >>> Extracting causal sites:
#> # of selected variants: 145

cat("Number of sites sampled with custom weights:", nrow(custom_sites$causal_genotypes), "\n")
#> Number of sites sampled with custom weights: 500
file.remove(tmp_custom_files)
#> [1] TRUE TRUE TRUE
```

------------------------------------------------------------------------

### Part 2: `simulate_phenotype()` - Specify the phenotype’s causal components

**Rationale:** Once a set of causal genotypes has been sampled, the
`simulate_phenotype()` function uses them to construct the final
phenotype. This function operationalizes the theoretical models
described earlier, allowing for the precise addition of genetic and
confounding effects.

#### Example 2.1: Basic quantitative trait

**Rationale:** The most fundamental simulation is of a trait with only a
genetic component ($g$) and a residual error component ($E$), with no
confounding. Notice that the population is highly structured, however,
there is no “environmental” confounding effect being simulated.
Therefore, the phenotypic covariance in this scenario can be fully
modeled via a “full GRM” that accounts for any level of relatedness as a
continuum (e.g. a `GeSi` matrix). This model serves as a baseline model
that let’s us to explicitly distingush the phenotypic covariance caused
by population structure, as opposed to environmental factors correlated
with genetic ancestry. Specifically, in this scenario, a full GRM is
enough to model the phenotypic covariance because there is no other
causal component. In practice, this means that a mixed linear model
(MLM) with a GRM but without principal components should be enough to
control for the inflation of the test statistics.

``` r
target_h2 <- 0.6

sim_basic <- simulate_phenotype(
  causal_genotypes = causal_geno_df, ## (p x n) data frame with genotypes at causal sites preceded by six metadata columns
  heritability = target_h2
)

# Validation: Check if the observed heritability is close to the target
observed_h2 <- var(sim_basic$phenotypes$g) / (var(sim_basic$phenotypes$g) + var(sim_basic$phenotypes$residual))
cat(sprintf("Target h2: %.3f | Observed h2: %.3f\n", target_h2, observed_h2))
#> Target h2: 0.600 | Observed h2: 0.596
```

#### Example 2.2: Adding a Categorical Confounder (Ordinal Scale)

**Rationale:** Phenocause can incorporate a categorical confounder
measured in either a nominal or ordinal scale. The scale is determined
by the class of the R vector. For instance, if the group assignment is
specified in the vector `categorical_confounder`, then a nominal scale
is assumed if the vector is of class `character`, and an ordinal scale
if it is of class `factor`. If the vector is of class `factor`, then the
confounding effects are monotonically increasing across categories in
the order given by `levels(categorical_confounder)`. Confounders in
ordinal scale can be useful to simulate variables with ordered levels,
like socioeconomic status or educational attainment.

In this example, we will create an ordinal confounder from the
population labels in the example metadata. However, notice that the
population levels themselves need not be the causal variables. Instead,
it there could be unknown causal variables differentially distributed
across the levels of whatever categorial variable was used to define the
groups (e.g. some unknown causal variables could be correlated with
educational attainment, which could in turn be correlated with ancestry
due to sociohistorical reasons).

We will order the population labels based on their mean value along the
first principal component (PC1). This simulates a situation where the
confounding effect is directionally associated with an axis of genetic
ancestry.

``` r
# Rationale: Create an ordered factor where the levels (population labels) are
# sorted by their average PC1 value. This establishes a monotonic relationship
# between the categorical variable and a major axis of genetic variation.
data(phenocause.metadata)
confounder_levels <- phenocause.metadata %>%
  group_by(population) %>%
  summarise(mean_pc1 = mean(PC1)) %>%
  arrange(mean_pc1) %>%
  pull(population)

categorical_confounder <- factor(phenocause.metadata$population, levels = confounder_levels)

## Note: if we wanted to use a nominal-scale confounder, we would simply do:
# categorical_confounder <- phenocause.metadata$population ## class "character"

cat("Ordered Levels for the Categorical Confounder:\n")
#> Ordered Levels for the Categorical Confounder:
print(levels(categorical_confounder))
#> [1] "MXB" "PEL" "CHB" "MXL" "IBS" "CLM" "PUR" "YRI"
```

Now, we simulate the phenotype, specifying the variance we want the
confounder to explain.

``` r
# Rationale: The `simulate_phenotype` function detects that the input is a
# factor and automatically simulates monotonic effects for its levels.
target_conf_var <- 0.1

sim_cat_conf <- simulate_phenotype(
  causal_genotypes = causal_geno_df,
  heritability = 0.5,
  categorical_confounder = categorical_confounder,
  categorical_confounder_variance = target_conf_var
)

# Validation: Check that the observed variance of the confounding effect
# is close to the target variance.
p_categorical <- sim_cat_conf$phenotypes$p_categorical_confounder
observed_conf_var <- var(p_categorical)

cat(sprintf("Target Confounder Variance: %.3f | Observed: %.3f\n", target_conf_var, observed_conf_var))
#> Target Confounder Variance: 0.100 | Observed: 0.100
```

We can also visualize the results to confirm that the effects are indeed
monotonic across the ordered levels.

``` r
# Rationale: Plotting the mean confounding effect for each category should
# reveal a clear, ordered trend, confirming the ordinal nature of the simulation.
plot_df <- data.frame(
  Confounder_Level = categorical_confounder,
  Confounding_Effect = p_categorical
)

ggplot(plot_df, aes(x = Confounder_Level, y = Confounding_Effect)) +
  geom_boxplot() +
  labs(
    title = "Ordinal Confounder Effects",
    x = "Confounder Level (Ordered by Mean PC1)",
    y = "Simulated Confounding Effect (P_cc)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![](man/figures/plot-ordinal-effects-1.png)<!-- -->

#### Example 2.3: Adding a confounder correlated with a specific ancestry component

**Rationale:** This simulates a common scenario in human genetics where
a non-genetic factor ($X_{ac}$) is correlated with genetic ancestry
($Q$) and also independently affects the phenotype. This allows for the
study of how historical and sociocultural variables that are correlated
with specific ancestry components can inflate GWAS test statistics and
bias downstream analyses. Here, we use the first principal component
(`PC1`) as our measure of ancestry, $Q$.

Notice that in this scenario the confounder effect is correlated with a
specific ancestry component as opposed to with the overall genetic
effect (See example 2.3). In this scenario, a GRM is not expected to be
enough to control the inflation of the test statistics in a MLM
association test.

``` r
set.seed(789)
# Load metadata containing sample info and PCs
data(phenocause.metadata)
pc1 <- phenocause.metadata$PC1

# Define simulation parameters
target_anc_cor <- 0.8 
target_anc_rel_var <- 0.15 

sim_anc_conf <- simulate_phenotype(
  causal_genotypes = causal_geno_df,
  heritability = 0.5,
  ancestry_component = pc1,
  ancestry_confounder_coeff = 0.6,
  ancestry_confounder_cor = target_anc_cor,
  ancestry_confounder_rel_var = target_anc_rel_var
)

# Validation: Check that the observed parameters match the targets
g_anc <- sim_anc_conf$phenotypes$g
p_ancestry <- sim_anc_conf$phenotypes$p_ancestry_confounder
x_ancestry <- sim_anc_conf$confounders$X_ancestry_confounder

observed_anc_cor <- cor(pc1, x_ancestry)
observed_anc_rel_var <- var(p_ancestry) / var(g_anc)

cat(sprintf("Target Cor(Q, X_ac): %.3f | Observed: %.3f\n", target_anc_cor, observed_anc_cor))
#> Target Cor(Q, X_ac): 0.800 | Observed: 0.799
cat(sprintf("Target Var(P_ac)/Var(g): %.3f | Observed: %.3f\n", target_anc_rel_var, observed_anc_rel_var))
#> Target Var(P_ac)/Var(g): 0.150 | Observed: 0.147
```

#### Example 2.4: Adding a confounder correlated with the polygenic score

**Rationale:** This models a distinct form of confounding where the
confounding variable ($X_{gc}$) is correlated with an individual’s true
polygenic score ($g$) itself.

``` r
set.seed(456)
target_g_cor <- 0.7 
target_g_rel_var <- 0.2

sim_genetic_conf <- simulate_phenotype(
  causal_genotypes = causal_geno_df,
  heritability = 0.5,
  genetic_confounder_coeff = 0.5,
  genetic_confounder_cor = target_g_cor,
  genetic_confounder_rel_var = target_g_rel_var
)

# Validation: Check the output
g <- sim_genetic_conf$phenotypes$g
p_genetic <- sim_genetic_conf$phenotypes$p_genetic_confounder
x_genetic <- sim_genetic_conf$confounders$X_genetic_confounder

observed_g_cor <- cor(g, x_genetic)
observed_g_rel_var <- var(p_genetic) / var(g)

cat(sprintf("Target Cor(g, X_gc): %.3f | Observed: %.3f\n", target_g_cor, observed_g_cor))
#> Target Cor(g, X_gc): 0.700 | Observed: 0.692
cat(sprintf("Target Var(P_gc)/Var(g): %.3f | Observed: %.3f\n", target_g_rel_var, observed_g_rel_var))
#> Target Var(P_gc)/Var(g): 0.200 | Observed: 0.201
```

#### Example 2.5: Full model for a binary trait & visualization

**Rationale:** We conclude by simulating a complex binary (case/control)
trait that incorporates all components: genetic effects and all three
types of confounding. For binary traits, we use the liability threshold
model. It is strongly recommended to set `liab_dist = "empirical"` when
including categorical confounders, as they can create a multimodal
liability distribution where the standard Gaussian assumption for the
threshold calculation would be inaccurate. We run the simulation in a
loop to assess the stability and robustness of the simulated parameters.

``` r
set.seed(1024)
# Define parameters 
target_prev <- 0.15
n_replicates <- 30 

# Define a categorical confounder (as a character vector for nominal effects)
categorical_confounder <- as.character(phenocause.metadata$population) 

# Run simulation across replicates
rep_results <- replicate(n_replicates, {
  sim_full <- simulate_phenotype(
    causal_genotypes =  causal_geno_df,
    heritability = 0.4,
    phenotype = "binary",
    prevalence = target_prev,
    liab_dist = "empirical",
    # Confounder 1: Categorical
    categorical_confounder = categorical_confounder,
    categorical_confounder_variance = 0.05,
    # Confounder 2: Ancestry-correlated
    ancestry_component = pc1,
    ancestry_confounder_coeff = 0.5,
    ancestry_confounder_rel_var = 0.1,
    ancestry_confounder_cor = 0.6,
    # Confounder 3: PGS-correlated
    genetic_confounder_coeff = -0.4,
    genetic_confounder_rel_var = 0.1,
    genetic_confounder_cor = -0.5
  )
  # Return key metrics for this replicate
  c(
    obs_prev = mean(sim_full$phenotypes$y_binary),
    obs_h2 = var(sim_full$phenotypes$g) / (var(sim_full$phenotypes$g) + var(sim_full$phenotypes$residual)),
    obs_anc_cor = cor(pc1, sim_full$confounders$X_ancestry_confounder),
    obs_anc_rel_var = var(sim_full$phenotypes$p_ancestry_confounder) / var(sim_full$phenotypes$g)
  )
}, simplify = FALSE)

results_df <- as.data.frame(do.call(rbind, rep_results))

# Summarize and check results
cat("Summary of metrics across", n_replicates, "replicates:\n")
#> Summary of metrics across 30 replicates:
print(summary(results_df))
#>     obs_prev        obs_h2        obs_anc_cor     obs_anc_rel_var  
#>  Min.   :0.15   Min.   :0.3935   Min.   :0.5850   Min.   :0.09591  
#>  1st Qu.:0.15   1st Qu.:0.3992   1st Qu.:0.5941   1st Qu.:0.09930  
#>  Median :0.15   Median :0.4021   Median :0.5992   Median :0.10022  
#>  Mean   :0.15   Mean   :0.4018   Mean   :0.5992   Mean   :0.10021  
#>  3rd Qu.:0.15   3rd Qu.:0.4053   3rd Qu.:0.6032   3rd Qu.:0.10092  
#>  Max.   :0.15   Max.   :0.4081   Max.   :0.6179   Max.   :0.10288

# Plot distributions of simulated parameter values vs. target values
p1 <- ggplot(results_df, aes(y = obs_prev)) +
  geom_boxplot() +
  geom_hline(yintercept = target_prev, color = "red", linetype = "dashed") +
  labs(title = "Prevalence", y = "Observed", x = "") +
  coord_flip() + theme_classic()

p2 <- ggplot(results_df, aes(y = obs_h2)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.4, color = "red", linetype = "dashed") +
  labs(title = "Heritability", y = "Observed", x = "") +
  coord_flip() + theme_classic()

p3 <- ggplot(results_df, aes(y = obs_anc_cor)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.6, color = "red", linetype = "dashed") +
  labs(title = "Cor(Q, Xac)", y = "Observed", x = "") +
  coord_flip() + theme_classic()

p4 <- ggplot(results_df, aes(y = obs_anc_rel_var)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed") +
  labs(title = "Var(Pac)/Var(g)", y = "Observed", x = "") +
  coord_flip() + theme_classic()

# Combine plots
p_merged <- (p1 | p2) / (p3 | p4) + 
  plot_annotation(title = paste("Simulation metrics across", n_replicates, "replicates"))

p_merged
```

![](man/figures/full-model-simulation-and-plot-1.png)<!-- -->

## References

1.  Baumdicker, F. et al. (2022). Efficient ancestry and mutation
    simulation with msprime 1.0. *Genetics*, 220, iyab229.  
2.  Chang, C. C. et al. (2015). Second-generation PLINK: rising to the
    challenge of larger and richer datasets. *GigaScience*, 4, 7.  
3.  Kelleher, J., Etheridge, A. M., & McVean, G. (2016). Efficient
    Coalescent Simulation and Genealogical Analysis for Large Sample
    Sizes. *PLOS Computational Biology*, 12, e1004842.  
4.  Medina-Muñoz, S. G. et al. (2023). Demographic modeling of admixed
    Latin American populations from whole genomes. *American Journal of
    Human Genetics*, 110(10), 1804–1816.  
5.  Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
    Improved Heritability Estimation from Genome-wide SNPs. *The
    American Journal of Human Genetics*, 91(6), 1011–1021.  
6.  Yang, J., Lee, S. H., Goddard, M. E., & Visscher, P. M. (2011).
    GCTA: A Tool for Genome-wide Complex Trait Analysis. *The American
    Journal of Human Genetics*, 88(1), 76–82.
