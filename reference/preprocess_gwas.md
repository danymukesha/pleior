# Preprocess GWAS Data

Filters and cleans GWAS summary statistics for downstream pleiotropy
analysis. Removes missing values, applies significance thresholds, and
standardizes column names.

## Usage

``` r
preprocess_gwas(
  gwas_data,
  pvalue_threshold = 5e-08,
  columns = c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG", "CHR_ID", "CHR_POS")
)
```

## Arguments

- gwas_data:

  A data.frame containing GWAS summary statistics.

- pvalue_threshold:

  Numeric. P-value threshold for filtering (default: 5e-8).

- columns:

  Character vector. Columns to retain in output (default: key GWAS
  columns).

## Value

A filtered and cleaned data.table.

## Examples

``` r
data(gwas_subset)
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
head(gwas_clean)
#>         SNPS          MAPPED_TRAIT PVALUE_MLOG CHR_ID   CHR_POS
#>       <char>                <char>       <num> <char>    <char>
#> 1:  rs814573     Alzheimer disease   672.69900     19  44908684
#> 2:  rs814573 myocardial infarction    15.00000     19  44908684
#> 3:    rs7412     Alzheimer disease   122.39790     19  44919689
#> 4:    rs7412       LDL cholesterol  9629.00000     19  44919689
#> 5: rs2817462    memory performance     6.30103      6 156588550
#> 6: rs2817462  tyrosine measurement     9.69897      6 156588550
```
