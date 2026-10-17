# GWAS Data for pleior Package

This directory contains real GWAS (Genome-Wide Association Study) data sourced from the official GWAS Catalog.

## Data Source

- **Source**: GWAS Catalog (European Bioinformatics Institute)
- **URL**: https://www.ebi.ac.uk/gwas/
- **Version**: 2024 release
- **License**: GWAS Catalog data is publicly available under CC0 1.0 Universal

## Data Files

### gwas_subset.rds
The example dataset in file form (same content as the `data(gwas_subset)`
package dataset):

- **Rows**: 2498 GWAS associations
- **Columns**: 5 essential columns for pleiotropy analysis
  - `SNPS`: SNP identifiers (rs numbers)
  - `MAPPED_TRAIT`: Standardized trait names
  - `PVALUE_MLOG`: -log10 transformed p-values
  - `CHR_ID`: Chromosome identifier (1-22, X, Y)
  - `CHR_POS`: Chromosomal position in base pairs

### example_gwas.tsv
A small tab-separated sample (first 200 rows) used to demonstrate
`load_gwas_data()`:

```r
file_path <- system.file("extdata", "example_gwas.tsv", package = "pleior")
gwas_data <- load_gwas_data(file_path)
```

### Traits Included

The dataset covers 15 traits:

1. Body mass index, body height
2. Alzheimer disease, schizophrenia
3. Myocardial infarction, coronary artery disease
4. Type 2 diabetes mellitus
5. Systolic / diastolic blood pressure
6. LDL, HDL and total cholesterol, triglycerides
7. C-reactive protein measurement
8. Asthma

### Data Processing Pipeline

The dataset was created using `data-raw/create_gwas_subset.R`:

1. **Load**: Read complete GWAS catalog associations file
2. **Select**: Keep associations for the 15 traits above
3. **Clean**: Standardize trait names and handle missing values
4. **Collapse**: One row per SNP-trait pair (most significant study kept)
5. **Sample**: Retain all associations of the strongest pleiotropic SNPs, then
   fill with a trait-balanced random sample to reach ~2500 rows
6. **Quality Check**: Ensure all required columns present and valid

## Pleiotropy in This Dataset

The curated subset contains **pleiotropic SNPs** - variants associated with multiple traits.

Example pleiotropic SNPs you may find:

- **rs429358**: Associated with 12 of the 15 traits (APOE region)
- **rs814573**: Associated with myocardial infarction and C-reactive protein
- **rs13107325**, **rs3184504**, **rs112875651** and others with multi-trait
  associations

## How to Use

Load the data in R:

```r
# Load pleior package
library(pleior)

# Load the example dataset
data(gwas_subset)

# View first few rows
head(gwas_subset)

# Basic statistics
cat("Total associations:", nrow(gwas_subset), "\n")
cat("Unique SNPs:", length(unique(gwas_subset$SNPS)), "\n")
cat("Unique traits:", length(unique(gwas_subset$MAPPED_TRAIT)), "\n")
```

## Complete Workflow Example

```r
library(pleior)

# Step 1: Load data
data(gwas_subset)

# Step 2: Preprocess - filter for genome-wide significance
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 5e-8)

# Step 3: Detect pleiotropic SNPs
pleio_results <- detect_pleiotropy(gwas_clean)

# Step 4: Explore results
cat("Found", nrow(pleio_results), "pleiotropic SNPs\n")
head(pleio_results)

# Step 5: Visualize
# Manhattan plot highlighting specific SNP
plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")

# Heatmap of SNP-trait associations
plot_pleiotropy_heatmap(pleio_results, top_n_snps = 20, top_n_traits = 10)

# Network visualization
plot_pleiotropy_network(pleio_results, top_n_snps = 15, show_labels = TRUE)
```

## Data Quality Notes

- All associations have significant p-values (PVALUE_MLOG >= 5, i.e., p < 1e-5)
- Chromosomal positions are validated
- Trait names are standardized from GWAS Catalog mapping
- SNP identifiers are cleaned (leading/trailing rs removed)
- Missing values handled appropriately

## Citation

When using this data in publications, please cite:

```
GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for genome-wide
association studies. *Nucleic Acids Research*, 51(D1), D944-D949.
https://doi.org/10.1093/nar/gkab1009
```

And cite the pleior package:

```
Mukesha, D. (2025). pleior: Pleiotropy Analysis for GWAS Data.
R package version 0.99.0. https://github.com/danymukesha/pleior
```

## Update Data

To update the dataset with newer GWAS catalog releases:

1. Download latest associations file from GWAS Catalog
2. Update `data-raw/gwas/gwas_associations.tsv`
3. Run `data-raw/create_real_dataset.R` to regenerate `data/gwas_subset.rda`
4. Update this README with any changes

## Additional Resources

- **GWAS Catalog**: https://www.ebi.ac.uk/gwas/
- **GWAS API**: https://www.ebi.ac.uk/gwas/api/docs/
- **Data Documentation**: https://www.ebi.ac.uk/gwas/docs/fileheaders/
