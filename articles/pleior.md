# Introduction to pleior

Abstract

Pleiotropy occurs when a single genetic variant influences multiple
traits, offering insights into shared biological pathways and disease
relationships. The package provides tools for analyzing pleiotropy in
genome-wide association studies (GWAS).

This vignette demonstrates a typical workflow using GWAS catalog data to
identify and visualize pleiotropic genetic variants.

## What is Pleiotropy?

Pleiotropy is a fundamental concept in genetics where a single gene or
genetic variant affects multiple phenotypic traits. In GWAS studies:

- **Genetic Pleiotropy**: One variant associated with multiple traits
- **Biological Significance**: Reveals shared pathways and mechanisms
- **Clinical Relevance**: Explains disease comorbidity and drug targets

Understanding pleiotropy helps the users to:

- Identify shared genetic architecture across diseases
- Discover novel therapeutic targets
- Explain why certain diseases co-occur
- Prioritize variants for functional studies

## Installation

``` r
# Install from GitHub
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("danymukesha/pleior")

# Load the package
library(pleior)
```

## Quick Start with Real GWAS Data

The package includes a curated subset of GWAS catalog data containing
~2,500 significant associations across 15 diverse traits.

``` r
# Load the example dataset
data(gwas_subset)

# Explore the data
head(gwas_subset)

# Summary statistics
cat("Dataset contains:", nrow(gwas_subset), "GWAS associations\n")
cat("Unique SNPs:", length(unique(gwas_subset$SNPS)), "\n")
cat("Unique traits:", length(unique(gwas_subset$MAPPED_TRAIT)), "\n")
```

## Complete Pleiotropy Analysis Workflow

### Step 1: Preprocess GWAS Data

Clean the data and filter for genome-wide significance:

``` r
# Filter for genome-wide significant associations (p < 5e-8)
gwas_clean <- preprocess_gwas(
    gwas_subset,
    pvalue_threshold = 5e-8
)

cat("After filtering:", nrow(gwas_clean), "significant associations\n")
```

### Step 2: Detect Pleiotropic SNPs

Identify SNPs associated with multiple traits:

``` r
# Detect pleiotropy across all traits
pleio_results <- detect_pleiotropy(gwas_clean)

cat("Found", nrow(pleio_results), "pleiotropic SNPs\n")

# Examine the top pleiotropic SNPs
head(pleio_results, 10)
```

### Step 3: Analyze Pleiotropic SNPs

``` r
# Summary statistics
snp_summary <- pleio_results %>%
    group_by(SNPS) %>%
    summarise(
        n_traits = n(),
        min_pvalue = max(PVALUE_MLOG),
        .groups = "drop"
    ) %>%
    arrange(desc(n_traits))

print(snp_summary)
```

### Step 4: Visualize Results

#### Manhattan Plot

Create a Manhattan plot highlighting pleiotropic SNPs:

``` r
# Plot with top pleiotropic SNP highlighted
if (nrow(pleio_results) > 0) {
    top_snp <- snp_summary$SNPS[1]

    p_manhattan <- plot_pleiotropy_manhattan(
        pleio_results,
        highlight_snp = top_snp,
        title = "Pleiotropic SNPs Across Genome"
    )

    print(p_manhattan)
}
```

#### Trait Distribution

See how many SNPs are associated with each trait:

``` r
trait_counts <- pleio_results %>%
    count(MAPPED_TRAIT, sort = TRUE)

# Simple bar plot
library(ggplot2)

p_traits <- ggplot(trait_counts, aes(x = reorder(MAPPED_TRAIT, n), y = n)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
        x = NULL,
        y = "Number of Pleiotropic SNPs",
        title = "Pleiotropic SNPs by Trait"
    ) +
    theme_minimal()

print(p_traits)
```

#### Pleiotropy Spectrum

Visualize the distribution of number of traits per SNP:

``` r
snp_trait_dist <- pleio_results %>%
    group_by(SNPS) %>%
    summarise(n_traits = n(), .groups = "drop")

p_spectrum <- ggplot(snp_trait_dist, aes(x = n_traits)) +
    geom_histogram(
        binwidth = 1,
        fill = "steelblue",
        color = "white",
        alpha = 0.8
    ) +
    labs(
        x = "Number of Traits per SNP",
        y = "Count",
        title = "Distribution of Pleiotropy"
    ) +
    theme_minimal()

print(p_spectrum)
```

## Advanced: Focus on Specific Traits

You can analyze pleiotropy between specific traits of interest:

``` r
# Define traits of interest
target_traits <- c(
    "Alzheimer's disease",
    "Type 2 diabetes",
    "Coronary artery disease"
)

# Filter for these traits
gwas_target <- gwas_clean %>%
    filter(MAPPED_TRAIT %in% target_traits)

# Detect pleiotropy among target traits
pleio_target <- detect_pleiotropy(
    gwas_target,
    traits = target_traits
)

cat("Pleiotropy between", length(target_traits), "target traits:\n")
print(pleio_target)
```

## Interpreting Results

### What Makes a SNP “Pleiotropic”?

A SNP is considered pleiotropic if it shows significant associations
with **two or more distinct traits** in our dataset. The package
identifies:

1.  **SNP Identifier**: The variant name (e.g., rs814573)
2.  **Associated Traits**: Which traits the SNP influences
3.  **Significance Level**: The strongest association (-log10 p-value)
4.  **Genomic Location**: Chromosome and position

### Example Interpretation

Consider a SNP like `rs814573`:

- **Associated with**: Alzheimer’s disease, myocardial infarction, LDL
  cholesterol
- **Significance**: p \< 5e-8 for all traits
- **Implication**: This variant affects pathways related to lipid
  metabolism and cardiovascular health
- **Research Direction**: Study shared biological mechanisms between
  these diseases

## Next Steps

After identifying pleiotropic SNPs:

1.  **Functional Annotation**: Map SNPs to genes and pathways
2.  **Literature Review**: Investigate known mechanisms
3.  **Experimental Validation**: Test in model systems
4.  **Drug Target Assessment**: Evaluate therapeutic potential

For advanced visualization methods including heatmaps, networks, and
publication-ready themes, see the **Visualization for Pleiotropy
Analysis** vignette.

## Performance Tips

For large GWAS datasets:

``` r
# Use data.table for efficient operations
library(data.table)

# Parallel processing can speed up analysis
library(future)
plan(multisession)

# Cache preprocessed results
if (file.exists("preprocessed_gwas.rds")) {
    gwas_clean <- readRDS("preprocessed_gwas.rds")
}
```

## Getting Help

``` r
# Function documentation
?detect_pleiotropy
?preprocess_gwas
?plot_pleiotropy_manhattan

# Package help
help(package = "pleior")
```

## References

- GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for
  genome-wide association studies. *Nucleic Acids Research*, 51(D1),
  D944-D949. <https://doi.org/10.1093/nar/gkab1009>

- Watanabe, K., et al. (2019). A global overview of pleiotropy and
  genetic architecture in complex traits. *Nature Genetics*, 51,
  1339-1348. <https://doi.org/10.1038/s41588-019-0685-1>

## Session Information

``` r
sessionInfo()
```
