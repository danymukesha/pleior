# Save Plot for Publication

Saves a ggplot2 object in publication-ready format with specified
resolution.

## Usage

``` r
save_publication_plot(
  plot,
  filename,
  width = 7,
  height = 5,
  dpi = 300,
  units = "in"
)
```

## Arguments

- plot:

  A ggplot2 object

- filename:

  Character. Output filename (extension determines format: .pdf, .png,
  .tiff, .svg)

- width:

  Numeric. Width in inches (default: 7)

- height:

  Numeric. Height in inches (default: 5)

- dpi:

  Numeric. Resolution for raster formats (default: 300 for publication)

- units:

  Character. Units: "in", "cm", "mm" (default: "in")

## Value

Invisible NULL (saves file to disk)

## Examples

``` r
if (FALSE) { # \dontrun{
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point()
save_publication_plot(p, "figure1.pdf", width = 7, height = 5)
} # }
```
