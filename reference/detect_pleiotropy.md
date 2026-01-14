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
#>          SNPS N_TRAITS                                  TRAITS
#>        <char>    <int>                                  <char>
#> 1: rs10182181        2    body mass index;tyrosine measurement
#> 2: rs10182181        2    body mass index;tyrosine measurement
#> 3:  rs2817462        2 memory performance;tyrosine measurement
#> 4:  rs3739081        2      C-reactive protein;body mass index
#> 5:  rs3739081        2      C-reactive protein;body mass index
#> 6:     rs7412        2       Alzheimer disease;LDL cholesterol
#>            MAPPED_TRAIT PVALUE_MLOG CHR_ID   CHR_POS
#>                  <char>       <num> <char>    <char>
#> 1:      body mass index    29.69897      2  24927427
#> 2: tyrosine measurement    11.69897      2  24927427
#> 3: tyrosine measurement     9.69897      6 156588550
#> 4:      body mass index     8.69897      2  26732753
#> 5:   C-reactive protein    11.00000      2  26732753
#> 6:    Alzheimer disease   122.39790     19  44919689

# Analyze specific traits
specific_traits <- c("Alzheimer disease", "myocardial infarction")
pleio_specific <- detect_pleiotropy(gwas_subset, traits = specific_traits)
```
