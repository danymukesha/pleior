# Detect Pleiotropic SNPs

Identifies SNPs associated with multiple traits in GWAS data. A SNP is
considered pleiotropic if it shows significant associations with two or
more distinct traits after applying the significance threshold.

## Usage

``` r
detect_pleiotropy(
  gwas_data,
  traits = NULL,
  pvalue_threshold = 5e-08,
  adjust_method = c("none", "bonferroni", "fdr"),
  clump_kb = NULL
)
```

## Arguments

- gwas_data:

  A data.frame containing GWAS summary statistics.

- traits:

  Character vector. Specific traits to analyze for pleiotropy (default:
  NULL uses all traits). Matching is exact against individual traits
  after splitting multi-trait fields.

- pvalue_threshold:

  Numeric. P-value threshold for significance (default: 5e-8).

- adjust_method:

  Character. Multiple-testing correction applied to all associations
  before thresholding: "none" (default), "bonferroni" or "fdr".

- clump_kb:

  Numeric or NULL. If supplied, co-located SNPs within this distance
  (kilobases) are clumped to their most significant representative. This
  is a position-based proxy for LD clumping (default: NULL, no
  clumping).

## Value

A data.frame with pleiotropic SNPs, one row per SNP-trait association,
including \`N_TRAITS\` and \`TRAITS\` summary columns.

## Details

The significance filter is applied \*before\* traits are counted, so
every reported pleiotropic SNP has at least two associations passing the
threshold. Multi-trait \`MAPPED_TRAIT\` fields (comma- or
semicolon-separated) are expanded into individual traits.

## Examples

``` r
data(gwas_subset)
pleio_results <- detect_pleiotropy(gwas_subset)
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

# Analyze specific traits
specific_traits <- c("Alzheimer disease", "myocardial infarction")
pleio_specific <- detect_pleiotropy(gwas_subset, traits = specific_traits)
```
