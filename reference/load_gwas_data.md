# Load GWAS Summary Statistics

Loads GWAS summary statistics from a file, supporting various formats.
The function automatically detects common separators and handles
standard GWAS file formats.

## Usage

``` r
load_gwas_data(file_path, sep = "\t", header = TRUE, quote = "")
```

## Arguments

- file_path:

  Character. Path to the GWAS summary statistics file.

- sep:

  Character. Field separator (default: "\t").

- header:

  Logical. Whether file contains header row (default: TRUE).

- quote:

  Character. Quote character (default: "").

## Value

A data.table containing the GWAS summary statistics.

## Examples

``` r
# \donttest{
# Load example data
file_path <- system.file("extdata", "example_gwas.tsv", package = "pleior")
if (file.exists(file_path)) {
    gwas_data <- load_gwas_data(file_path)
    head(gwas_data)
}
# }
```
