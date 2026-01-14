# Trait Co-occurrence Network Visualization

Creates a network showing how often pairs of traits share pleiotropic
SNPs. This helps identify trait clusters and shared biological
mechanisms.

## Usage

``` r
plot_trait_cooccurrence_network(
  pleio_data,
  min_shared_snps = 1,
  top_n_traits = 15,
  node_size = 8,
  layout = "fr",
  show_edge_labels = TRUE,
  color_by = "n_snps"
)
```

## Arguments

- pleio_data:

  A data.frame containing pleiotropy analysis results

- min_shared_snps:

  Integer. Minimum shared SNPs to show edge (default: 1)

- top_n_traits:

  Integer. Number of top traits to include (default: 15)

- node_size:

  Numeric. Base node size (default: 8)

- layout:

  Character. Layout algorithm: "fr", "kk", "dh", "grid", "circle"
  (default: "fr")

- show_edge_labels:

  Logical. Show shared SNP counts on edges (default: TRUE)

- color_by:

  Character. Color nodes by: "n_snps", "degree", "none" (default:
  "n_snps")

## Value

A ggraph ggplot2 object

## Examples

``` r
data(gwas_subset)
pleio_results <- detect_pleiotropy(gwas_subset)
if (nrow(pleio_results) > 0) {
    p <- plot_trait_cooccurrence_network(pleio_results, top_n_traits = 5)
    print(p)
}
#> Warning: No trait pairs meet minimum shared SNP criteria

```
