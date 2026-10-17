# Adjust p-values for Multiple Testing

Converts the stored -log10(p-values), applies a multiple-testing
correction across all supplied tests, and returns the -log10 transformed
adjusted p-values in the same column.

## Usage

``` r
adjust_gwas_pvalues(gwas_data, method = "none", pvalue_col = "PVALUE_MLOG")
```

## Arguments

- gwas_data:

  A data.frame.

- method:

  Character. One of "none", "bonferroni" or "fdr".

- pvalue_col:

  Character. Name of the -log10(p) column (default: "PVALUE_MLOG").

## Value

A data.frame with adjusted -log10(p-values).
