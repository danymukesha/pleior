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
#>           SNPS                                    MAPPED_TRAIT PVALUE_MLOG
#>         <char>                                          <char>       <num>
#> 1: rs115866895                                     body height   13.698970
#> 2: rs560221037 low density lipoprotein cholesterol measurement   11.221849
#> 3:  rs12081179                         systolic blood pressure    7.522879
#> 4:   rs4908677                                     body height    7.522879
#> 5:   rs6682378                                 body mass index    9.698970
#> 6:   rs1381928                                          asthma    8.397940
#>    CHR_ID CHR_POS
#>     <int>   <int>
#> 1:      1 1661177
#> 2:      1 2235046
#> 3:      1 2319680
#> 4:      1 7678120
#> 5:      1 8071697
#> 6:      1 8630874
# }
```
