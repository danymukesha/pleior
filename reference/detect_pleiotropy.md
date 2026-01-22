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
#>   SNPS       N_TRAITS TRAITS   `DISEASE/TRAIT` CHR_ID CHR_POS `REPORTED GENE(S)`
#>   <chr>         <int> <chr>    <chr>           <chr>  <chr>   <chr>             
#> 1 rs10006235        3 Alzheim… Educational at… 4      129748… Intergenic        
#> 2 rs10006766        2 hemoglo… Protein quanti… 4      872388… NR                
#> 3 rs10007186        2 creatin… Non-albumin pr… 4      786678… ANXA3             
#> 4 rs1000778        25 HbA1c m… Sphingolipid l… 11     618878… FADS3             
#> 5 rs10008637       13 cholest… Estimated glom… 4      764929… SHROOM3           
#> 6 rs10008637       13 cholest… Cardiometaboli… 4      764929… SHROOM3           
#> # ℹ 4 more variables: `P-VALUE` <chr>, PVALUE_MLOG <dbl>,
#> #   `P-VALUE (TEXT)` <chr>, MAPPED_TRAIT <chr>

# Analyze specific traitsdetect_pleiotropydetect_pleiotropy
specific_traits <- c("Alzheimer disease", "myocardial infarction")
pleio_specific <- detect_pleiotropy(gwas_subset, traits = specific_traits)
```
