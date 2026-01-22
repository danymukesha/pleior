# Regional Association Plot (Zoomed-in View)

Creates a regional association plot showing detailed view of genomic
region around a specific pleiotropic SNP. Shows all associations in the
region across all traits, highlighting the pleiotropic SNP of interest.

## Usage

``` r
plot_regional_association(
  pleio_data,
  target_snp,
  window_size = 5e+05,
  highlight_color = "#E31A1C",
  show_gene_names = FALSE,
  gene_data = NULL
)
```

## Arguments

- pleio_data:

  A data.frame containing pleiotropy analysis results

- target_snp:

  Character. Target SNP ID to zoom in on (required)

- window_size:

  Integer. Window size in base pairs (default: 500000)

- highlight_color:

  Character. Color for target SNP (default: "#E31A1C")

- show_gene_names:

  Logical. Show nearby gene names (default: FALSE)

- gene_data:

  Optional data.frame with gene annotations (CHR, START, END, NAME)

## Value

A ggplot2 regional plot object

## Examples

``` r
data(gwas_subset)
pleio_results <- detect_pleiotropy(gwas_subset)
if (nrow(pleio_results) > 0 && "rs814573" %in% pleio_results$SNPS) {
    p <- plot_regional_association(pleio_results, target_snp = "rs814573")
    print(p)
}
```
