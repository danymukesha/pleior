# Split Multi-Trait Fields into One Row per Trait

GWAS Catalog \`MAPPED_TRAIT\` values may contain several traits
separated by commas or semicolons (e.g. "body mass index, obesity").
Counting such a field as a single trait understates pleiotropy, so this
helper explodes it into one row per individual trait.

## Usage

``` r
split_mapped_traits(gwas_data, trait_col = "MAPPED_TRAIT")
```

## Arguments

- gwas_data:

  A data.frame.

- trait_col:

  Character. Name of the trait column (default: "MAPPED_TRAIT").

## Value

A data.frame with one row per (association, individual trait).
