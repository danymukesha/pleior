# Detect Pleiotropic SNPs

Identifies SNPs associated with multiple traits in GWAS data. A SNP is
considered pleiotropic if it shows significant associations with two or
more distinct traits.

## Usage

``` r
detect_pleiotropy(gwas_data, traits = NULL, pvalue_threshold = 5e-08)
```

## Arguments

- gwas_data:

  A data.frame containing preprocessed GWAS summary statistics.

- traits:

  Character vector. Specific traits to analyze for pleiotropy (default:
  NULL uses all traits).

- pvalue_threshold:

  Numeric. P-value threshold for significance (default: 5e-8).

## Value

A data.table with pleiotropic SNPs, their associated traits, and
significance levels.

## Examples

``` r
data(gwas_subset)
pleio_results <- detect_pleiotropy(gwas_subset)
head(pleio_results)
#> # A tibble: 6 × 11
#>   SNPS         N_TRAITS TRAITS `DISEASE/TRAIT` CHR_ID CHR_POS `REPORTED GENE(S)`
#>   <chr>           <int> <chr>  <chr>           <chr>  <chr>   <chr>             
#> 1 1:221765779…        2 force… Lung function … NA     NA      NA                
#> 2 1:237929787…        2 force… Lung function … NA     NA      RYR2              
#> 3 6p21.32             2 disor… Chronic inflam… NA     NA      NA                
#> 4 6p21.32             2 disor… Pharyngeal dis… NA     NA      NA                
#> 5 6p21.32             2 disor… Sinonasal dise… NA     NA      NA                
#> 6 A*01:01             3 susce… Shingles        NA     NA      HLA-A             
#> # ℹ 4 more variables: `P-VALUE` <chr>, PVALUE_MLOG <dbl>,
#> #   `P-VALUE (TEXT)` <chr>, MAPPED_TRAIT <chr>

# Analyze specific traits
specific_traits <- c("Alzheimer disease", "myocardial infarction")
pleio_specific <- detect_pleiotropy(gwas_subset, traits = specific_traits)
```
