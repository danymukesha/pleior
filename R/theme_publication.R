#' Publication-Ready Theme for Pleiotropy Plots
#'
#' Creates a publication-ready ggplot2 theme optimized for scientific journals.
#' Supports multiple journal styles and includes colorblind-friendly palettes.
#'
#' @param journal_style Character. Journal style to use. Options: "nature", "science", "pnas", "default" (default: "default")
#' @param base_size Numeric. Base font size (default: 12)
#' @param legend_position Character. Legend position: "right", "bottom", "none" (default: "right")
#'
#' @return A ggplot2 theme object
#'
#' @import ggplot2
#' @importFrom scales viridis_pal
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point()
#' p + theme_pleiotropy_publication(journal_style = "nature")
#'
#' @export
theme_pleiotropy_publication <- function(journal_style = "default",
                                         base_size = 12,
                                         legend_position = "right") {
    theme_base <- theme_minimal(base_size = base_size) +
        theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
            panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.text = element_text(color = "black", size = base_size * 0.8),
            axis.title = element_text(color = "black", size = base_size * 0.9, face = "bold"),
            plot.title = element_text(color = "black", size = base_size * 1.1, face = "bold", hjust = 0.5),
            legend.position = legend_position,
            legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.25),
            legend.title = element_text(face = "bold", size = base_size * 0.8),
            legend.text = element_text(size = base_size * 0.75),
            strip.background = element_rect(fill = "gray90", color = "gray80"),
            strip.text = element_text(face = "bold", size = base_size * 0.8)
        )

    journal_specific <- switch(journal_style,
        "nature" = theme(
            axis.title = element_text(size = 8, face = "plain"),
            axis.text = element_text(size = 7),
            plot.title = element_text(size = 9, face = "plain", hjust = 0),
            legend.position = "none",
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 6),
            aspect.ratio = 0.8
        ),
        "science" = theme(
            axis.title = element_text(size = 10, face = "bold"),
            axis.text = element_text(size = 9),
            plot.title = element_text(size = 12, face = "bold", hjust = 0),
            legend.position = "right",
            aspect.ratio = 0.75
        ),
        "pnas" = theme(
            axis.title = element_text(size = 9, face = "bold"),
            axis.text = element_text(size = 8),
            plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            legend.position = "right",
            legend.title = element_text(size = 8),
            aspect.ratio = 0.85
        ),
        theme()
    )

    theme_base + journal_specific
}

#' Colorblind-Friendly Palettes for Pleiotropy Visualization
#'
#' Returns colorblind-friendly color palettes for visualization.
#'
#' @param palette_name Character. Palette name: "okabe_ito", "viridis", "plasma", "cividis", "blue_red", "default" (default: "okabe_ito")
#' @param n Integer. Number of colors to return (default: NULL returns all colors)
#'
#' @return Character vector of hex color codes
#'
#' @importFrom scales viridis_pal
#'
#' @examples
#' get_pleiotropy_colors("okabe_ito", n = 5)
#' get_pleiotropy_colors("viridis", n = 10)
#'
#' @export
get_pleiotropy_colors <- function(palette_name = "okabe_ito", n = NULL) {
    palettes <- list(
        "okabe_ito" = c(
            "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#000000", "#999999", "#F0F0F0"
        ),
        "viridis" = scales::viridis_pal(option = "D")(10),
        "plasma" = scales::viridis_pal(option = "C")(10),
        "cividis" = scales::viridis_pal(option = "E")(10),
        "blue_red" = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"),
        "default" = c(
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
        )
    )

    if (!palette_name %in% names(palettes)) {
        stop("Invalid palette_name. Choose from: ", paste(names(palettes), collapse = ", "))
    }

    colors <- palettes[[palette_name]]

    if (!is.null(n)) {
        colors <- colorRampPalette(colors)(n)
    }

    return(colors)
}

#' Save Plot for Publication
#'
#' Saves a ggplot2 object in publication-ready format with specified resolution.
#'
#' @param plot A ggplot2 object
#' @param filename Character. Output filename (extension determines format: .pdf, .png, .tiff, .svg)
#' @param width Numeric. Width in inches (default: 7)
#' @param height Numeric. Height in inches (default: 5)
#' @param dpi Numeric. Resolution for raster formats (default: 300 for publication)
#' @param units Character. Units: "in", "cm", "mm" (default: "in")
#'
#' @return Invisible NULL (saves file to disk)
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point()
#' save_publication_plot(p, "figure1.pdf", width = 7, height = 5)
#' }
#'
#' @export
save_publication_plot <- function(plot, filename, width = 7, height = 5,
                                  dpi = 300, units = "in") {
    if (!inherits(plot, "ggplot")) {
        stop("plot must be a ggplot2 object")
    }

    ext <- tolower(tools::file_ext(filename))

    if (!ext %in% c("pdf", "png", "tiff", "svg", "eps")) {
        stop("Unsupported file format. Use: .pdf, .png, .tiff, .svg, or .eps")
    }

    if (ext == "pdf") {
        ggsave(filename,
            plot = plot, width = width, height = height,
            units = units, device = "pdf"
        )
    } else if (ext == "png") {
        ggsave(filename,
            plot = plot, width = width, height = height,
            units = units, dpi = dpi
        )
    } else if (ext == "tiff") {
        ggsave(filename,
            plot = plot, width = width, height = height,
            units = units, dpi = dpi
        )
    } else if (ext == "svg") {
        ggsave(filename,
            plot = plot, width = width, height = height,
            units = units, device = "svg"
        )
    } else if (ext == "eps") {
        ggsave(filename,
            plot = plot, width = width, height = height,
            units = units, device = "eps"
        )
    }

    invisible(NULL)
}

#' Combine Multiple Plots for Publication
#'
#' Combines multiple ggplot2 objects into a single figure using patchwork.
#' Supports grid layouts and labels for multi-panel figures.
#'
#' @param ... Multiple ggplot2 objects
#' @param ncol Integer. Number of columns (default: NULL for auto)
#' @param nrow Integer. Number of rows (default: NULL for auto)
#' @param labels Character vector. Panel labels (e.g., c("A", "B", "C"))
#' @param tag_position Character. Position of panel labels: "topleft", "top", "topright", etc. (default: "topleft")
#'
#' @return A patchwork object
#'
#' @importFrom patchwork wrap_plots plot_layout
#'
#' @examples
#' \dontrun{
#' p1 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point()
#' p2 <- ggplot(mtcars, aes(x = wt, y = hp)) +
#'     geom_point()
#' combine_publication_plots(p1, p2, ncol = 2, labels = c("A", "B"))
#' }
#'
#' @export
combine_publication_plots <- function(..., ncol = NULL, nrow = NULL,
                                      labels = NULL, tag_position = "topleft") {
    plots <- list(...)

    if (length(plots) == 0) {
        stop("No plots provided")
    }

    combined <- wrap_plots(plots, ncol = ncol, nrow = nrow)

    if (!is.null(labels) && length(labels) == length(plots)) {
        combined <- combined + plot_layout(tag_level = tag_position) &
            theme(plot.tag = element_text(size = 12, face = "bold"))

        for (i in seq_along(plots)) {
            combined <- combined + patchwork::patchworkGrob(
                ggplot2::annotation_custom(
                    grid::textGrob(labels[i],
                        x = grid::unit(0.05, "npc"),
                        y = grid::unit(0.95, "npc"), gp = grid::gpar(fontsize = 12, fontface = "bold")
                    )
                )
            )
        }
    }

    return(combined)
}
