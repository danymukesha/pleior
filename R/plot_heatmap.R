#' Heatmap Visualization for SNP-Trait Associations
#'
#' Creates a publication-ready heatmap showing -log10(p-value) associations between
#' pleiotropic SNPs and traits. Useful for visualizing patterns of shared genetics.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param top_n_snps Integer. Number of top pleiotropic SNPs to include (default: 20)
#' @param top_n_traits Integer. Number of top traits to include (default: 10)
#' @param value_col Character. Column to use for heatmap values (default: "PVALUE_MLOG")
#' @param show_dendrograms Logical. Show row and column dendrograms (default: FALSE)
#' @param clustering_method Character. Clustering method: "complete", "average", "ward.D2" (default: "complete")
#' @param color_palette Character. Color palette: "bluered", "viridis", "plasma", "default" (default: "bluered")
#' @param scale Character. Scale values: "none", "row", "column" (default: "none")
#' @param show_values Logical. Show numeric values in cells (default: FALSE)
#' @param value_format Character. Format string for values (default: "%.1f")
#' @param legend_title Character. Legend title (default: "-log10(p)")
#'
#' @return A ggplot2 heatmap object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter arrange slice_head ungroup select rename
#' @importFrom tidyr pivot_wider
#' @importFrom scales viridis_pal
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_pleiotropy_heatmap(pleio_results, top_n_snps = 5, top_n_traits = 4)
#'     print(p)
#' }
#'
#' @export
plot_pleiotropy_heatmap <- function(pleio_data,
                                    top_n_snps = 20,
                                    top_n_traits = 10,
                                    value_col = "PVALUE_MLOG",
                                    show_dendrograms = FALSE,
                                    clustering_method = "complete",
                                    color_palette = "bluered",
                                    scale = "none",
                                    show_values = FALSE,
                                    value_format = "%.1f",
                                    legend_title = "-log10(p)") {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (!value_col %in% names(pleio_data)) {
        stop("Column '", value_col, "' not found in data")
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    snp_counts <- pleio_data |>
        group_by(SNPS) |>
        summarise(N_TRAITS = n(), .groups = "drop") |>
        arrange(desc(N_TRAITS)) |>
        slice_head(n = top_n_snps)

    trait_counts <- pleio_data |>
        group_by(MAPPED_TRAIT) |>
        summarise(N_SNPS = n(), .groups = "drop") |>
        arrange(desc(N_SNPS)) |>
        slice_head(n = top_n_traits)

    plot_data <- pleio_data |>
        filter(SNPS %in% snp_counts$SNPS, MAPPED_TRAIT %in% trait_counts$MAPPED_TRAIT) |>
        group_by(SNPS, MAPPED_TRAIT) |>
        summarise(value = max(!!sym(value_col)), .groups = "drop") |>
        ungroup()

    if (nrow(plot_data) == 0) {
        warning("No data after filtering")
        return(ggplot() +
            theme_void() +
            annotate("text",
                x = 0.5, y = 0.5,
                label = "No data to display\nAdjust filtering criteria"
            ))
    }

    heatmap_data <- plot_data |>
        pivot_wider(names_from = MAPPED_TRAIT, values_from = value, values_fill = 0) |>
        mutate(across(where(is.numeric), ~ ifelse(is.na(.x), 0, .x)))

    snp_order <- heatmap_data$SNPS
    trait_order <- colnames(heatmap_data)[-1]

    if (!show_dendrograms) {
        snp_counts_ordered <- plot_data |>
            group_by(SNPS) |>
            summarise(total = sum(value), .groups = "drop") |>
            arrange(desc(total))

        trait_counts_ordered <- plot_data |>
            group_by(MAPPED_TRAIT) |>
            summarise(total = sum(value), .groups = "drop") |>
            arrange(desc(total))

        heatmap_data$SNPS <- factor(heatmap_data$SNPS, levels = snp_counts_ordered$SNPS)

        plot_data_long <- heatmap_data |>
            pivot_longer(cols = -SNPS, names_to = "MAPPED_TRAIT", values_to = "value") |>
            mutate(MAPPED_TRAIT = factor(MAPPED_TRAIT, levels = trait_counts_ordered$MAPPED_TRAIT))
    } else {
        snp_dist <- dist(heatmap_data[, -1], method = "euclidean")
        snp_hc <- hclust(snp_dist, method = clustering_method)

        trait_dist <- dist(t(heatmap_data[, -1]), method = "euclidean")
        trait_hc <- hclust(trait_dist, method = clustering_method)

        heatmap_data$SNPS <- factor(heatmap_data$SNPS, levels = rownames(heatmap_data)[snp_hc$order])

        plot_data_long <- heatmap_data |>
            pivot_longer(cols = -SNPS, names_to = "MAPPED_TRAIT", values_to = "value") |>
            mutate(MAPPED_TRAIT = factor(MAPPED_TRAIT, levels = colnames(heatmap_data)[-1][trait_hc$order]))
    }

    colors <- switch(color_palette,
        "bluered" = c("#0571b0", "#92c5de", "#f7f7f7", "#f4a582", "#ca0020"),
        "viridis" = scales::viridis_pal(option = "D")(100),
        "plasma" = scales::viridis_pal(option = "C")(100),
        "default" = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"),
        scales::viridis_pal(option = "D")(100)
    )

    p <- ggplot(plot_data_long, aes(x = MAPPED_TRAIT, y = SNPS, fill = value)) +
        geom_tile(color = "white", size = 0.5) +
        scale_fill_gradientn(
            colors = colors,
            name = legend_title,
            guide = guide_colorbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = unit(15, "lines"),
                barheight = unit(0.5, "lines")
            )
        ) +
        labs(
            title = "Pleiotropic SNP-Trait Association Heatmap",
            subtitle = paste(
                "Top", nlevels(plot_data_long$SNPS), "SNPs ×",
                nlevels(plot_data_long$MAPPED_TRAIT), "Traits"
            ),
            x = "Trait",
            y = "SNP"
        ) +
        coord_fixed() +
        theme_pleiotropy_publication() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            legend.position = "bottom"
        )

    if (show_values) {
        p <- p + geom_text(aes(label = sprintf(value_format, value)),
            size = 2, color = "black"
        )
    }

    return(p)
}

#' Regional Association Plot (Zoomed-in View)
#'
#' Creates a regional association plot showing detailed view of genomic region
#' around a specific pleiotropic SNP. Shows all associations in the region across
#' all traits, highlighting the pleiotropic SNP of interest.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param target_snp Character. Target SNP ID to zoom in on (required)
#' @param window_size Integer. Window size in base pairs (default: 500000)
#' @param highlight_color Character. Color for target SNP (default: "#E31A1C")
#' @param show_gene_names Logical. Show nearby gene names (default: FALSE)
#' @param gene_data Optional data.frame with gene annotations (CHR, START, END, NAME)
#'
#' @return A ggplot2 regional plot object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter ungroup select
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0 && "rs814573" %in% pleio_results$SNPS) {
#'     p <- plot_regional_association(pleio_results, target_snp = "rs814573")
#'     print(p)
#' }
#'
#' @export
plot_regional_association <- function(pleio_data,
                                      target_snp,
                                      window_size = 500000,
                                      highlight_color = "#E31A1C",
                                      show_gene_names = FALSE,
                                      gene_data = NULL) {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (missing(target_snp)) {
        stop("target_snp is required")
    }

    required_cols <- c("SNPS", "CHR_ID", "CHR_POS", "PVALUE_MLOG", "MAPPED_TRAIT")
    missing_cols <- setdiff(required_cols, names(pleio_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    target_data <- pleio_data |>
        filter(SNPS == target_snp)

    if (nrow(target_data) == 0) {
        stop("Target SNP not found in data")
    }

    target_chr <- target_data$CHR_ID[1]
    target_pos <- as.numeric(target_data$CHR_POS[1])

    region_start <- max(0, target_pos - window_size)
    region_end <- target_pos + window_size

    region_data <- pleio_data |>
        filter(
            CHR_ID == target_chr,
            as.numeric(CHR_POS) >= region_start,
            as.numeric(CHR_POS) <= region_end
        ) |>
        mutate(
            CHR_POS = as.numeric(CHR_POS),
            distance_kb = (CHR_POS - target_pos) / 1000,
            is_target = SNPS == target_snp
        ) |>
        ungroup()

    if (nrow(region_data) == 0) {
        warning("No SNPs in the specified window")
        return(ggplot() +
            theme_void())
    }

    colors <- get_pleiotropy_colors("viridis", n = length(unique(region_data$MAPPED_TRAIT)))

    p <- ggplot(region_data, aes(x = distance_kb, y = PVALUE_MLOG, color = MAPPED_TRAIT)) +
        geom_point(alpha = 0.7, size = 1.5) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        scale_color_manual(values = colors, name = "Trait") +
        labs(
            title = paste("Regional Association Plot:", target_snp),
            subtitle = sprintf("Chr %s: %d ± %d kb", target_chr, target_pos, window_size / 1000),
            x = "Distance from target SNP (kb)",
            y = expression(-log[10](p))
        ) +
        theme_pleiotropy_publication() +
        theme(
            legend.position = "right",
            panel.grid.minor = element_blank()
        )

    target_snp_data <- region_data |>
        filter(is_target)

    if (nrow(target_snp_data) > 0) {
        p <- p +
            geom_point(
                data = target_snp_data,
                color = highlight_color,
                size = 3,
                shape = 17,
                alpha = 0.9
            ) +
            annotate("text",
                x = 0, y = max(region_data$PVALUE_MLOG) * 0.95,
                label = target_snp, color = highlight_color,
                size = 4, fontface = "bold", hjust = 0.5
            )
    }

    if (show_gene_names && !is.null(gene_data)) {
        gene_regions <- gene_data |>
            filter(
                CHR == target_chr,
                START >= region_start,
                END <= region_end
            ) |>
            mutate(
                mid_pos = (START + END) / 2,
                distance_kb = (mid_pos - target_pos) / 1000
            )

        if (nrow(gene_regions) > 0) {
            p <- p +
                annotate("rect",
                    xmin = (gene_regions$START - target_pos) / 1000,
                    xmax = (gene_regions$END - target_pos) / 1000,
                    ymin = -Inf, ymax = -Inf,
                    alpha = 0.2, fill = "gray70"
                ) +
                annotate("text",
                    x = gene_regions$distance_kb,
                    y = min(region_data$PVALUE_MLOG) * 0.05,
                    label = gene_regions$NAME,
                    angle = 45, hjust = 0, size = 2.5, color = "gray30"
                )
        }
    }

    return(p)
}
