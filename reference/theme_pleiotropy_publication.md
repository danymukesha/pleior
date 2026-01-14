# Publication-Ready Theme for Pleiotropy Plots

Creates a publication-ready ggplot2 theme optimized for scientific
journals. Supports multiple journal styles and includes
colorblind-friendly palettes.

## Usage

``` r
theme_pleiotropy_publication(
  journal_style = "default",
  base_size = 12,
  legend_position = "right"
)
```

## Arguments

- journal_style:

  Character. Journal style to use. Options: "nature", "science", "pnas",
  "default" (default: "default")

- base_size:

  Numeric. Base font size (default: 12)

- legend_position:

  Character. Legend position: "right", "bottom", "none" (default:
  "right")

## Value

A ggplot2 theme object

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point()
p + theme_pleiotropy_publication(journal_style = "nature")

```
