# Combine Multiple Plots for Publication

Combines multiple ggplot2 objects into a single figure using patchwork.
Supports grid layouts and labels for multi-panel figures.

## Usage

``` r
combine_publication_plots(
  ...,
  ncol = NULL,
  nrow = NULL,
  labels = NULL,
  tag_position = "topleft"
)
```

## Arguments

- ...:

  Multiple ggplot2 objects

- ncol:

  Integer. Number of columns (default: NULL for auto)

- nrow:

  Integer. Number of rows (default: NULL for auto)

- labels:

  Character vector. Panel labels (e.g., c("A", "B", "C"))

- tag_position:

  Character. Position of panel labels: "topleft", "top", "topright",
  etc. (default: "topleft")

## Value

A patchwork object

## Examples

``` r
if (FALSE) { # \dontrun{
p1 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
    geom_point()
p2 <- ggplot(mtcars, aes(x = wt, y = hp)) +
    geom_point()
combine_publication_plots(p1, p2, ncol = 2, labels = c("A", "B"))
} # }
```
