#' Venn Diagram for SNP Sharing Between Traits
#'
#' Creates Venn diagrams showing overlap of pleiotropic SNPs between 2-4 traits.
#' Helps visualize shared genetic architecture across phenotypes.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param traits Character vector. 2-4 traits to compare (required)
#' @param title Character. Plot title (default: "SNP Sharing Between Traits")
#' @param show_counts Logical. Show counts in each region (default: TRUE)
#' @param show_percentages Logical. Show percentages (default: FALSE)
#' @param color_palette Character. Color palette: "okabe_ito", "viridis", "default" (default: "okabe_ito")
#' @param alpha Numeric. Transparency (default: 0.5)
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @importFrom dplyr filter group_by summarise mutate
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_venn_diagram(pleio_results,
#'         traits = c("Alzheimer disease", "myocardial infarction")
#'     )
#'     print(p)
#' }
#'
#' @export
plot_venn_diagram <- function(pleio_data,
                              traits,
                              title = "SNP Sharing Between Traits",
                              show_counts = TRUE,
                              show_percentages = FALSE,
                              color_palette = "okabe_ito",
                              alpha = 0.5) {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (missing(traits) || length(traits) < 2 || length(traits) > 4) {
        stop("Provide 2-4 traits to compare")
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    trait_snps <- lapply(traits, function(t) {
        unique(pleio_data$SNPS[pleio_data$MAPPED_TRAIT == t])
    })
    names(trait_snps) <- traits

    colors <- get_pleiotropy_colors(color_palette, n = length(traits))

    if (length(traits) == 2) {
        venn_plot <- plot_venn_2set(
            trait_snps, colors, alpha,
            show_counts, show_percentages, title
        )
    } else if (length(traits) == 3) {
        venn_plot <- plot_venn_3set(
            trait_snps, colors, alpha,
            show_counts, show_percentages, title
        )
    } else if (length(traits) == 4) {
        venn_plot <- plot_venn_4set(
            trait_snps, colors, alpha,
            show_counts, show_percentages, title
        )
    }

    return(venn_plot)
}

plot_venn_2set <- function(trait_snps, colors, alpha,
                           show_counts, show_percentages, title) {
    snp1 <- trait_snps[[1]]
    snp2 <- trait_snps[[2]]

    only1 <- length(setdiff(snp1, snp2))
    only2 <- length(setdiff(snp2, snp1))
    intersection <- length(intersect(snp1, snp2))
    total <- length(unique(c(snp1, snp2)))

    label1 <- if (show_counts) {
        paste0(only1, if (show_percentages) paste0("\n(", round(only1 / total * 100, 1), "%)", ""))
    } else if (show_percentages) {
        paste0(round(only1 / total * 100, 1), "%")
    } else {
        ""
    }

    label2 <- if (show_counts) {
        paste0(only2, if (show_percentages) paste0("\n(", round(only2 / total * 100, 1), "%)", ""))
    } else if (show_percentages) {
        paste0(round(only2 / total * 100, 1), "%")
    } else {
        ""
    }

    label_intersect <- if (show_counts) {
        paste0(intersection, if (show_percentages) paste0("\n(", round(intersection / total * 100, 1), "%)", ""))
    } else if (show_percentages) {
        paste0(round(intersection / total * 100, 1), "%")
    } else {
        ""
    }

    df <- data.frame(
        x = c(0.35, 0.65),
        y = c(0.5, 0.5),
        r = 0.35,
        label = c(label1, label2),
        color = colors[1:2],
        name = names(trait_snps)
    )

    p <- ggplot() +
        geom_circle(aes(x0 = x, y0 = y, r = r, fill = color, color = color),
            data = df, alpha = alpha, linewidth = 1
        ) +
        geom_text(aes(x = x, y = y + 0.35, label = name),
            data = df,
            size = 4, fontface = "bold"
        ) +
        geom_text(aes(x = 0.45, y = 0.5, label = label_intersect), size = 5, fontface = "bold") +
        scale_fill_identity() +
        scale_color_identity() +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        labs(title = title) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
        annotate("text", x = 0.2, y = 0.5, label = label1, size = 5, fontface = "bold") +
        annotate("text", x = 0.8, y = 0.5, label = label2, size = 5, fontface = "bold")

    return(p)
}

plot_venn_3set <- function(trait_snps, colors, alpha,
                           show_counts, show_percentages, title) {
    snp1 <- trait_snps[[1]]
    snp2 <- trait_snps[[2]]
    snp3 <- trait_snps[[3]]

    total <- length(unique(c(snp1, snp2, snp3)))

    only1 <- length(setdiff(snp1, union(snp2, snp3)))
    only2 <- length(setdiff(snp2, union(snp1, snp3)))
    only3 <- length(setdiff(snp3, union(snp1, snp2)))
    intersect12 <- length(intersect(setdiff(snp1, snp3), setdiff(snp2, snp3)))
    intersect13 <- length(intersect(setdiff(snp1, snp2), setdiff(snp3, snp2)))
    intersect23 <- length(intersect(setdiff(snp2, snp1), setdiff(snp3, snp1)))
    intersect_all <- length(Reduce(intersect, list(snp1, snp2, snp3)))

    df <- data.frame(
        x = c(0.5, 0.35, 0.65),
        y = c(0.55, 0.3, 0.3),
        r = 0.3,
        color = colors[1:3],
        name = names(trait_snps)
    )

    positions <- list(
        only1 = c(0.5, 0.65),
        only2 = c(0.3, 0.3),
        only3 = c(0.7, 0.3),
        intersect12 = c(0.4, 0.45),
        intersect13 = c(0.6, 0.45),
        intersect23 = c(0.5, 0.3),
        intersect_all = c(0.5, 0.45)
    )

    counts <- c(only1, only2, only3, intersect12, intersect13, intersect23, intersect_all)

    labels <- sapply(counts, function(n) {
        if (show_counts) {
            paste0(n, if (show_percentages) paste0("\n(", round(n / total * 100, 1), "%)", ""))
        } else if (show_percentages) {
            paste0(round(n / total * 100, 1), "%")
        } else {
            ""
        }
    })

    p <- ggplot() +
        geom_circle(aes(x0 = x, y0 = y, r = r, fill = color, color = color),
            data = df, alpha = alpha, linewidth = 1
        ) +
        geom_text(aes(x = x, y = y + 0.35, label = name),
            data = df,
            size = 3.5, fontface = "bold"
        ) +
        mapply(function(pos, lab) {
            annotate("text", x = pos[1], y = pos[2], label = lab, size = 4, fontface = "bold")
        }, positions, labels) +
        scale_fill_identity() +
        scale_color_identity() +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        labs(title = title) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    return(p)
}

plot_venn_4set <- function(trait_snps, colors, alpha,
                           show_counts, show_percentages, title) {
    snp1 <- trait_snps[[1]]
    snp2 <- trait_snps[[2]]
    snp3 <- trait_snps[[3]]
    snp4 <- trait_snps[[4]]

    total <- length(unique(c(snp1, snp2, snp3, snp4)))

    only1 <- length(setdiff(snp1, Reduce(union, list(snp2, snp3, snp4))))
    only2 <- length(setdiff(snp2, Reduce(union, list(snp1, snp3, snp4))))
    only3 <- length(setdiff(snp3, Reduce(union, list(snp1, snp2, snp4))))
    only4 <- length(setdiff(snp4, Reduce(union, list(snp1, snp2, snp3))))

    intersect_all <- length(Reduce(intersect, list(snp1, snp2, snp3, snp4)))

    labels <- c(only1, only2, only3, only4, rep(intersect_all, 11))
    positions <- list(
        only1 = c(0.25, 0.65),
        only2 = c(0.75, 0.65),
        only3 = c(0.25, 0.25),
        only4 = c(0.75, 0.25),
        intersect12 = c(0.5, 0.65),
        intersect13 = c(0.25, 0.45),
        intersect14 = c(0.4, 0.45),
        intersect23 = c(0.75, 0.45),
        intersect24 = c(0.6, 0.45),
        intersect34 = c(0.5, 0.25),
        intersect123 = c(0.35, 0.55),
        intersect124 = c(0.65, 0.55),
        intersect134 = c(0.35, 0.35),
        intersect234 = c(0.65, 0.35),
        intersect_all = c(0.5, 0.45)
    )

    df <- data.frame(
        x = c(0.3, 0.7, 0.3, 0.7),
        y = c(0.6, 0.6, 0.3, 0.3),
        r = 0.25,
        color = colors[1:4],
        name = names(trait_snps)
    )

    label_text <- sapply(labels[1:15], function(n) {
        if (show_counts) {
            paste0(n, if (show_percentages) paste0("\n(", round(n / total * 100, 1), "%)", ""))
        } else if (show_percentages) {
            paste0(round(n / total * 100, 1), "%")
        } else {
            ""
        }
    })

    p <- ggplot() +
        geom_circle(aes(x0 = x, y0 = y, r = r, fill = color, color = color),
            data = df, alpha = alpha, linewidth = 1
        ) +
        geom_text(aes(x = x, y = y + 0.3, label = name),
            data = df,
            size = 3, fontface = "bold"
        ) +
        mapply(function(pos, lab) {
            annotate("text", x = pos[1], y = pos[2], label = lab, size = 3, fontface = "bold")
        }, positions[1:15], label_text) +
        scale_fill_identity() +
        scale_color_identity() +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        labs(title = title) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    return(p)
}

#' Pleiotropy Significance Landscape
#'
#' Creates a cumulative distribution plot showing the significance landscape
#' of pleiotropic associations. Helps understand the distribution of
#' significance across all pleiotropic SNPs.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param n_bins Integer. Number of bins for histogram (default: 30)
#' @param show_median Logical. Show median line (default: TRUE)
#' @param show_mean Logical. Show mean line (default: FALSE)
#' @param show_quantiles Logical. Show 25th and 75th percentiles (default: TRUE)
#' @param color_palette Character. Color palette (default: "viridis")
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_pleiotropy_landscape(pleio_results)
#'     print(p)
#' }
#'
#' @export
plot_pleiotropy_landscape <- function(pleio_data,
                                      n_bins = 30,
                                      show_median = TRUE,
                                      show_mean = FALSE,
                                      show_quantiles = TRUE,
                                      color_palette = "viridis") {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (!"PVALUE_MLOG" %in% names(pleio_data)) {
        stop("Column 'PVALUE_MLOG' not found in data")
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    pvalues <- pleio_data$PVALUE_MLOG
    median_val <- median(pvalues, na.rm = TRUE)
    mean_val <- mean(pvalues, na.rm = TRUE)
    q25 <- quantile(pvalues, 0.25, na.rm = TRUE)
    q75 <- quantile(pvalues, 0.75, na.rm = TRUE)

    colors <- get_pleiotropy_colors(color_palette, n = 10)

    p <- ggplot(data.frame(PVALUE_MLOG = pvalues), aes(x = PVALUE_MLOG)) +
        geom_histogram(aes(y = after_stat(density)),
            bins = n_bins,
            fill = colors[1],
            color = colors[1],
            alpha = 0.6
        ) +
        geom_density(aes(color = "Density"), linewidth = 1.5, alpha = 0) +
        scale_color_manual(values = "black", name = "") +
        labs(
            title = "Pleiotropy Significance Landscape",
            subtitle = sprintf("N = %d associations", length(pvalues)),
            x = expression(-log[10](p)),
            y = "Density"
        ) +
        theme_pleiotropy_publication() +
        theme(
            legend.position = "top",
            panel.grid.minor = element_blank()
        )

    if (show_median) {
        p <- p +
            geom_vline(
                xintercept = median_val,
                linetype = "dashed",
                color = colors[3],
                linewidth = 1,
                alpha = 0.8
            ) +
            annotate("text",
                x = median_val, y = Inf,
                label = sprintf("Median: %.2f", median_val),
                vjust = 2, hjust = -0.05, size = 3.5, color = colors[3]
            )
    }

    if (show_mean) {
        p <- p +
            geom_vline(
                xintercept = mean_val,
                linetype = "dotdash",
                color = colors[2],
                linewidth = 1,
                alpha = 0.8
            ) +
            annotate("text",
                x = mean_val, y = Inf,
                label = sprintf("Mean: %.2f", mean_val),
                vjust = 3.5, hjust = -0.05, size = 3.5, color = colors[2]
            )
    }

    if (show_quantiles) {
        p <- p +
            geom_vline(
                xintercept = q25,
                linetype = "dotted",
                color = "gray50",
                linewidth = 0.8
            ) +
            geom_vline(
                xintercept = q75,
                linetype = "dotted",
                color = "gray50",
                linewidth = 0.8
            ) +
            annotate("rect",
                xmin = q25, xmax = q75, ymin = -Inf, ymax = Inf,
                alpha = 0.1, fill = colors[1]
            )
    }

    return(p)
}
