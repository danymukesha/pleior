# GWAS Summary Statistics Subset

A curated subset of GWAS Catalog summary statistics prepared for
pleiotropy analysis.

## Usage

``` r
gwas_subset
```

## Format

A data.table with 2498 rows and 5 columns:

- SNPS:

  Character. SNP identifiers (rs numbers)

- MAPPED_TRAIT:

  Character. Standardized trait or phenotype name

- PVALUE_MLOG:

  Numeric. -log10 transformed p-values for significance

- CHR_ID:

  Character. Chromosome identifier (1-22, X, Y)

- CHR_POS:

  Numeric. Chromosomal position in base pairs

## Source

GWAS Catalog associations file (2024 release), downloaded from:
<https://www.ebi.ac.uk/gwas/api/search/downloads/alternative>

The subset was created by `data-raw/create_gwas_subset.R`: 1. Loading
the complete GWAS Catalog associations file 2. Selecting 15 traits
relevant to complex disease and metabolism 3. Collapsing duplicate study
entries to one row per SNP-trait pair (keeping the most significant
association) 4. Retaining all associations of the strongest pleiotropic
SNPs, then filling with a trait-balanced random sample to reach ~2500
rows 5. Ensuring representation across multiple chromosomes

## Details

This dataset contains genome-wide association study (GWAS) summary
statistics from the GWAS Catalog (https://www.ebi.ac.uk/gwas/), reduced
to a small, reproducible set that demonstrates pleiotropy detection and
visualization. It is the example dataset shipped with \`pleior\`.

The dataset contains 2498 associations covering 1835 distinct SNPs
across 15 traits. Of these SNPs, 128 are pleiotropic (associated with
two or more of the 15 traits); the strongest, `rs429358` (in the APOE
region), maps to 12 traits. `rs814573`, used as a running example in the
README, is associated with myocardial infarction and C-reactive protein
measurement in this subset.

The 15 traits are:

- body mass index

- body height

- Alzheimer disease

- myocardial infarction

- coronary artery disease

- type 2 diabetes mellitus

- systolic blood pressure

- diastolic blood pressure

- low density lipoprotein cholesterol measurement

- high density lipoprotein cholesterol measurement

- triglyceride measurement

- total cholesterol measurement

- C-reactive protein measurement

- asthma

- schizophrenia

The data is suitable for:

- Demonstrating pleiotropy detection (SNPs associated with multiple
  traits)

- Testing visualization functions

- Teaching GWAS analysis workflows

- Benchmarking analysis methods

## References

GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for genome-wide
association studies. \*Nucleic Acids Research\*, 51(D1), D944-D949.
<https://doi.org/10.1093/nar/gkab1009>

## Examples

``` r
# Load the data
data(gwas_subset)
head(gwas_subset)
#>           SNPS                                    MAPPED_TRAIT PVALUE_MLOG
#>         <char>                                          <char>       <num>
#> 1: rs115866895                                     body height   13.698970
#> 2: rs560221037 low density lipoprotein cholesterol measurement   11.221849
#> 3:  rs12081179                         systolic blood pressure    7.522879
#> 4:   rs4908677                                     body height    7.522879
#> 5:   rs6682378                                 body mass index    9.698970
#> 6:   rs1381928                                          asthma    8.397940
#>    CHR_ID CHR_POS
#>    <char>   <num>
#> 1:      1 1661177
#> 2:      1 2235046
#> 3:      1 2319680
#> 4:      1 7678120
#> 5:      1 8071697
#> 6:      1 8630874

# Basic statistics
n_unique_snps <- length(unique(gwas_subset$SNPS))
n_unique_traits <- length(unique(gwas_subset$MAPPED_TRAIT))
cat("Unique SNPs:", n_unique_snps, "\n")
#> Unique SNPs: 1835 
cat("Unique traits:", n_unique_traits, "\n")
#> Unique traits: 15 

# Preprocess and detect pleiotropy
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
pleio_results <- detect_pleiotropy(gwas_clean)
head(pleio_results)
#>       SNPS                                     MAPPED_TRAIT PVALUE_MLOG CHR_ID
#> 1 rs429358                   C-reactive protein measurement  1952.69897     19
#> 2 rs429358 high density lipoprotein cholesterol measurement   497.00000     19
#> 3 rs429358                                Alzheimer disease   323.00000     19
#> 4 rs429358                    total cholesterol measurement   323.00000     19
#> 5 rs429358  low density lipoprotein cholesterol measurement   186.52288     19
#> 6 rs429358                          coronary artery disease    42.39794     19
#>    CHR_POS N_TRAITS
#> 1 44908684       12
#> 2 44908684       12
#> 3 44908684       12
#> 4 44908684       12
#> 5 44908684       12
#> 6 44908684       12
#>                                                                                                                                                                                                                                                                                                                                TRAITS
#> 1 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
#> 2 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
#> 3 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
#> 4 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
#> 5 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
#> 6 Alzheimer disease;C-reactive protein measurement;body height;body mass index;coronary artery disease;high density lipoprotein cholesterol measurement;low density lipoprotein cholesterol measurement;myocardial infarction;systolic blood pressure;total cholesterol measurement;triglyceride measurement;type 2 diabetes mellitus
```
