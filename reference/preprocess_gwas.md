# Preprocess GWAS Data

Filters and cleans GWAS summary statistics for downstream pleiotropy
analysis. Removes missing values, applies significance thresholds,
optionally adjusts for multiple testing, and standardizes column names.

## Usage

``` r
preprocess_gwas(
  gwas_data,
  pvalue_threshold = 5e-08,
  columns = c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG", "CHR_ID", "CHR_POS"),
  adjust_method = c("none", "bonferroni", "fdr")
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

- adjust_method:

  Character. Multiple-testing correction applied before thresholding:
  "none" (default), "bonferroni" or "fdr".

## Value

A filtered and cleaned data.frame (of the same class as the input).

## Details

The -log10(p-value) column (\`PVALUE_MLOG\`) is derived automatically
from a raw p-value column (\`P-VALUE\`, \`PVALUE\` or \`P\`) when it is
not present.

## Examples

``` r
data(gwas_subset)
gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
head(gwas_clean)
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
```
