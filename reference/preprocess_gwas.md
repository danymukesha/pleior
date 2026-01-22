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
#>          SNPS    MAPPED_TRAIT PVALUE_MLOG CHR_ID   CHR_POS
#>        <char>          <char>       <num> <char>    <char>
#> 1:  rs2820292 body mass index    10.09691      1 201815159
#> 2:  rs6725549 body mass index    74.00000      2    630323
#> 3: rs10182181 body mass index    29.69897      2  24927427
#> 4:  rs3739081 body mass index     8.69897      2  26732753
#> 5: rs13011109 body mass index    14.00000      2  58630284
#> 6:  rs1016287 body mass index    12.39794      2  59078490
```
