# Colorblind-Friendly Palettes for Pleiotropy Visualization

Returns colorblind-friendly color palettes for visualization.

## Usage

``` r
get_pleiotropy_colors(palette_name = "okabe_ito", n = NULL)
```

## Arguments

- palette_name:

  Character. Palette name: "okabe_ito", "viridis", "plasma", "cividis",
  "blue_red", "default" (default: "okabe_ito")

- n:

  Integer. Number of colors to return (default: NULL returns all colors)

## Value

Character vector of hex color codes

## Examples

``` r
get_pleiotropy_colors("okabe_ito", n = 5)
#> [1] "#E69F00" "#3CAF66" "#6A6859" "#321E29" "#F0F0F0"
get_pleiotropy_colors("viridis", n = 10)
#>  [1] "#440154" "#482878" "#3E4A89" "#31688E" "#26828E" "#1F9E89" "#35B779"
#>  [8] "#6DCD59" "#B4DE2C" "#FDE725"
```
