# pleior <a href="https://danymukesha.github.io/pleior/"><img src="man/figures/logo.png" align="right" height="139" alt="pleior website" /></a>

This package is designed to analyze pleiotropy in genome-wide association
studies (GWAS). Pleiotropy occurs when a single genetic variant influences
multiple traits, offering insights into shared biological pathways and disease
relationships. This package enables users to identify, analyze, and visualize
pleiotropic associations, addressing question of how genetic variants
contribute to multiple complex traits, such as Alzheimer's disease and related
conditions.

## Why It Matters

Understanding pleiotropy is crucial for uncovering shared biological
mechanisms, explaining disease comorbidity, and identifying potential
therapeutic targets. However, existing tools for pleiotropy analysis are often
fragmented, focusing on specific aspects like detection or visualization, and
may lack user-friendly interfaces or integration with functional annotations.

Recent research suggests that pleiotropy is widespread in human genetics, with
many genetic variants linked to multiple diseases or traits. For example, SNP
(rs814573) is known to be associated with Alzheimer's disease and other traits
like myocardial infarction and tyrosine measurement. Understanding these
connections could help explain why certain diseases co-occur (co-morbidity) to
identify potential targets and responses of treatment/drug based on shared
genetic predispositions with other known traits. `pleior` provides a unified
platform to explore these relationships, making it easier to generate
hypotheses and interpret results.

A 2019 study in *Nature Genetics* found that trait-associated loci cover over
half of genome, with 90% overlapping multiple traits. Similarly, studies like
FactorGo (2023) and PLACO (2020) emphasize the need for scalable,
statistically robust methods to characterize pleiotropy across thousands of
traits. `pleior` builds on these findings by providing a unified toolset to
explore pleiotropy in a practical, accessible manner.

## Installation

### Development Version

To install the development version from GitHub:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("danymukesha/pleior")
```

## Quick Start

```r
library(pleior)

# Load example dataset (real GWAS catalog data)
data(gwas_subset)

# View dataset summary
cat("Dataset contains:", nrow(gwas_subset), "GWAS associations\n")
cat("Unique SNPs:", length(unique(gwas_subset$SNPS)), "\n")
cat("Unique traits:", length(unique(gwas_subset$MAPPED_TRAIT)), "\n")
head(gwas_subset)
```

## GWAS Data

The package includes a **curated subset of GWAS Catalog data** containing
approximately 2,500 significant associations across 15 diverse traits.

### Data Source

- **Source**: GWAS Catalog (European Bioinformatics Institute)
- **Version**: 2024 release
- **URL**: https://www.ebi.ac.uk/gwas/
- **License**: GWAS Catalog data is publicly available under CC0 1.0 Universal

### Included Traits

The dataset includes diverse traits to demonstrate pleiotropy:

- Body mass index (BMI)
- Alzheimer's disease
- Type 2 diabetes
- Coronary artery disease
- Blood pressure traits
- Lipid levels (LDL, HDL cholesterol)
- Other metabolic traits

### Data Quality

- All associations have significant p-values (p < 1e-5)
- Chromosomal positions validated
- Trait names standardized from GWAS Catalog mapping
- SNP identifiers cleaned and validated

## Features

### Data Loading and Preprocessing

- Load GWAS summary statistics from various formats
- Filter and standardize GWAS data
- Handle missing values and data validation
- Support for genome-wide significance filtering

### Pleiotropy Detection

- Identify SNPs associated with multiple traits
- Analyze specific trait combinations
- Calculate association strength
- Generate pleiotropy statistics

### Visualization

- **Manhattan plots**: Genome-wide association visualization
- **Heatmaps**: SNP-trait association matrices
- **Network graphs**: Trait-SNP relationship networks
- **Venn diagrams**: Set visualization for trait overlaps
- **Regional association plots**: Zoomed views around specific SNPs
- **Publication-ready themes**: Optimized for major journals
- **Colorblind-friendly palettes**: Okabe-Ito and other accessible palettes

## Complete Workflow Example

```r
library(pleior)

# Step 1: Load GWAS data
data(gwas_subset)

# Step 2: Preprocess - filter for genome-wide significance
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 5e-8)

# Step 3: Detect pleiotropic SNPs
pleio_results <- detect_pleiotropy(gwas_clean)

# Step 4: Explore results
cat("Found", nrow(pleio_results), "pleiotropic SNPs\n")
head(pleio_results)

# Step 5: Visualize results

# Manhattan plot highlighting specific SNP
plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")

# Heatmap of SNP-trait associations
plot_pleiotropy_heatmap(pleio_results, top_n_snps = 20, top_n_traits = 10)

# Network visualization
plot_pleiotropy_network(pleio_results, top_n_snps = 15, show_labels = TRUE)

# Venn diagram comparing traits
traits_to_compare <- unique(pleio_results$MAPPED_TRAIT)[1:3]
plot_venn_diagram(pleio_results, traits = traits_to_compare)
```

## Understanding the Results

### What Makes a SNP Pleiotropic?

In this dataset, a SNP is considered pleiotropic if it shows significant
associations with **two or more distinct traits**. The output includes:

- **SNP Identifier**: The variant name (e.g., rs814573)
- **Associated Traits**: Which traits the SNP influences
- **Significance Level**: The strongest association (-log10 p-value)
- **Genomic Location**: Chromosome and position

### Example: rs814573

This SNP is associated with:

- Alzheimer's disease
- Myocardial infarction
- LDL cholesterol

**Implications**:

- Affects pathways related to lipid metabolism and cardiovascular health
- Potential drug target for multiple conditions
- Reveals shared genetic mechanisms between diseases

## Documentation

Detailed documentation is available through:

### Vignettes

- **Introduction to pleior**: Complete workflow with GWAS data
- **Visualization for Pleiotropy Analysis**: Advanced visualization techniques

Access vignettes:

```r
vignette("Introduction to pleior", package = "pleior")
vignette("Visualization for Pleiotropy Analysis", package = "pleior")
```

### Function Documentation

Each function has detailed help pages:

```r
# Core functions
?detect_pleiotropy
?preprocess_gwas
?load_gwas_data

# Visualization functions
?plot_pleiotropy_manhattan
?plot_pleiotropy_heatmap
?plot_pleiotropy_network
?plot_venn_diagram

# Utility functions
?theme_pleiotropy_publication
?get_pleiotropy_colors
```

## Citation

If you use `pleior` in published research, please cite:

### Package Citation
```
Mukesha, D. (2025). pleior: Pleiotropy Analysis for GWAS Data.
R package version 0.99.0. https://github.com/danymukesha/pleior
```

Generate citation within R:
```r
citation("pleior")
```

### Data Citation

The included GWAS dataset should be cited as:
```
GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for genome-wide
association studies. *Nucleic Acids Research*, 51(D1), D944-D949.
https://doi.org/10.1093/nar/gkab1009
```

## Performance Tips

For analyzing large GWAS datasets:

1. **Use data.table**: Efficient handling of large files
2. **Filter early**: Reduce dataset size early in workflow
3. **Parallel processing**: Use future package for parallelization
4. **Cache results**: Save intermediate results to avoid recomputation

## Troubleshooting

### Common Issues

**Issue**: "No pleiotropic SNPs found"

- **Solution**: Lower the p-value threshold in `preprocess_gwas()` or check
data quality

**Issue**: "Plot takes too long to render"

- **Solution**: Reduce top_n parameters or subset data before plotting

**Issue**: Memory errors with large datasets

- **Solution**: Process data in chunks or use filtering to reduce size

## Support

- **Issues**: Report bugs at https://github.com/danymukesha/pleior/issues
- **Documentation**: See package vignettes and function help pages
- **GWAS Catalog**: https://www.ebi.ac.uk/gwas/docs/

## Acknowledgments

This package uses data and resources from:

- GWAS Catalog (European Bioinformatics Institute)
- GWAS Catalog API
- R community packages (ggplot2, dplyr, igraph, etc.)
