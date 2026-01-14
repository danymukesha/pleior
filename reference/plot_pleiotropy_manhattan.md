# Create Manhattan Plot for Pleiotropic SNPs

Generates a Manhattan plot highlighting pleiotropic SNPs across
chromosomes. The plot shows the distribution of significant associations
and can highlight specific SNPs of interest.

## Usage

``` r
plot_pleiotropy_manhattan(
  pleio_data,
  highlight_snp = NULL,
  title = "Manhattan Plot of Pleiotropic SNPs"
)
```

## Arguments

- pleio_data:

  A data.frame containing pleiotropy analysis results.

- highlight_snp:

  Character. SNP identifier to highlight (default: NULL).

- title:

  Character. Plot title (default: "Manhattan Plot of Pleiotropic SNPs").

## Value

A ggplot2 object representing the Manhattan plot.

## Examples

``` r
data(gwas_subset)
pleio_results <- detect_pleiotropy(gwas_subset)
if (nrow(pleio_results) > 0) {
    p <- plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")
    print(p)
}

```
