#' Multi-Trait Manhattan Plot
#'
#' Creates a faceted Manhattan plot showing genome-wide associations for multiple
#' traits simultaneously. Each panel represents a different trait, making it easy
#' to compare patterns across traits.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param traits Character vector. Specific traits to plot (default: NULL uses top traits)
#' @param max_traits Integer. Maximum number of traits to display (default: 6)
#' @param significance_line Numeric. Genome-wide significance threshold (default: 5e-8)
#' * @param suggestive_line Numeric. Suggestive significance threshold (default: 1e-5)
#' @param highlight_snp Character. SNP to highlight across all traits (default: NULL)
#' @param point_size Numeric. Point size (default: 1)
#' @param alpha Numeric. Point transparency (default: 0.6)
#' @param ncol Integer. Number of columns in facet grid (default: 2)
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter arrange slice_head ungroup left_join
#' @importFrom tidyr unite
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_multi_trait_manhattan(pleio_results, max_traits = 3)
#'     print(p)
#' }
#'
#' @export
plot_multi_trait_manhattan <- function(pleio_data,
                                       traits = NULL,
                                       max_traits = 6,
                                       significance_line = 5e-8,
                                       suggestive_line = 1e-5,
                                       highlight_snp = NULL,
                                       point_size = 1,
                                       alpha = 0.6,
                                       ncol = 2) {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    required_cols <- c("SNPS", "CHR_ID", "CHR_POS", "PVALUE_MLOG", "MAPPED_TRAIT")
    missing_cols <- setdiff(required_cols, names(pleio_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    if (is.null(traits)) {
        trait_counts <- pleio_data |>
            group_by(MAPPED_TRAIT) |>
            summarise(N_SNPS = n(), .groups = "drop") |>
            arrange(desc(N_SNPS)) |>
            slice_head(n = max_traits)

        traits <- trait_counts$MAPPED_TRAIT
    }

    plot_data <- pleio_data |>
        filter(MAPPED_TRAIT %in% traits) |>
        mutate(
            CHR_ID = as.character(CHR_ID),
            CHR_POS = as.numeric(CHR_POS),
            CHR_NUM = as.numeric(ifelse(CHR_ID == "X", "23",
                ifelse(CHR_ID == "Y", "24", CHR_ID)
            ))
        ) |>
        filter(!is.na(CHR_NUM), !is.na(CHR_POS)) |>
        arrange(MAPPED_TRAIT, CHR_NUM, CHR_POS)

    if (nrow(plot_data) == 0) {
        warning("No data after filtering")
        return(ggplot() +
            theme_void())
    }

    chr_lengths <- plot_data |>
        group_by(MAPPED_TRAIT, CHR_ID, CHR_NUM) |>
        summarise(MAX_POS = max(CHR_POS), .groups = "drop") |>
        arrange(MAPPED_TRAIT, CHR_NUM) |>
        group_by(MAPPED_TRAIT) |>
        mutate(
            CUMPOS_START = cumsum(c(0, MAX_POS[-n()])),
            CUMPOS_END = cumsum(MAX_POS),
            CENTER = CUMPOS_START + (MAX_POS / 2)
        ) |>
        ungroup()

    plot_data <- plot_data |>
        left_join(chr_lengths[, c("MAPPED_TRAIT", "CHR_ID", "CUMPOS_START")],
            by = c("MAPPED_TRAIT", "CHR_ID")
        ) |>
        mutate(PLOT_POS = CUMPOS_START + CHR_POS)

    colors <- get_pleiotropy_colors("okabe_ito", n = 24)

    sig_threshold <- -log10(significance_line)
    sugg_threshold <- -log10(suggestive_line)

    p <- ggplot(plot_data, aes(
        x = PLOT_POS, y = PVALUE_MLOG,
        color = CHR_ID
    )) +
        geom_point(size = point_size, alpha = alpha) +
        geom_hline(
            yintercept = sig_threshold,
            linetype = "solid", color = "red",
            linewidth = 0.5, alpha = 0.7
        ) +
        geom_hline(
            yintercept = sugg_threshold,
            linetype = "dashed", color = "blue",
            linewidth = 0.5, alpha = 0.7
        ) +
        facet_wrap(~MAPPED_TRAIT, ncol = ncol, scales = "free_x") +
        scale_color_manual(values = colors) +
        labs(
            title = "Multi-Trait Manhattan Plot",
            subtitle = paste0(
                "Significance threshold: ", significance_line, " (red), ",
                suggestive_line, " (blue)"
            ),
            x = "Chromosome",
            y = expression(-log[10](p))
        ) +
        theme_pleiotropy_publication() +
        theme(
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            strip.text = element_text(face = "bold", size = 10)
        )

    facet_data <- chr_lengths |>
        group_by(MAPPED_TRAIT) |>
        summarise(
            chr_labels = paste(CHR_ID, collapse = ";"),
            chr_centers = paste(CENTER, collapse = ";"),
            .groups = "drop"
        )

    plot_data <- plot_data |>
        left_join(facet_data, by = "MAPPED_TRAIT")

    p <- p + scale_x_continuous(
        expand = c(0.01, 0),
        labels = waiver(),
        breaks = waiver()
    )

    if (!is.null(highlight_snp)) {
        highlight_data <- plot_data |>
            filter(SNPS == highlight_snp)

        if (nrow(highlight_data) > 0) {
            p <- p +
                geom_point(
                    data = highlight_data,
                    color = "darkred",
                    size = point_size * 2,
                    alpha = 1, shape = 17
                ) +
                geom_text(
                    data = highlight_data,
                    aes(label = SNPS),
                    vjust = -0.5, color = "darkred", size = 2.5,
                    check_overlap = TRUE
                )
        }
    }

    return(p)
}

#' Effect Size Comparison Plot
#'
#' Creates a comparison plot showing effect sizes across traits for pleiotropic SNPs.
#' Helps identify which SNPs have stronger effects on specific traits.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results
#' @param top_n_snps Integer. Number of top pleiotropic SNPs to include (default: 15)
#' @param effect_col Character. Column to use for effect size (default: "PVALUE_MLOG")
#' @param use_log_scale Logical. Use log scale for effect sizes (default: TRUE)
#' @param show_error_bars Logical. Show error bars (default: FALSE)
#' @param error_col Character. Column for standard errors (default: NULL)
#' @param color_by Character. Color by: "snp", "trait", "chromosome" (default: "snp")
#' @param flip_coords Logical. Flip x and y coordinates (default: TRUE for horizontal bars)
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate filter arrange slice_head ungroup
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_effect_size_comparison(pleio_results, top_n_snps = 5)
#'     print(p)
#' }
#'
#' @export
plot_effect_size_comparison <- function(pleio_data,
                                        top_n_snps = 15,
                                        effect_col = "PVALUE_MLOG",
                                        use_log_scale = TRUE,
                                        show_error_bars = FALSE,
                                        error_col = NULL,
                                        color_by = "snp",
                                        flip_coords = TRUE) {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (!effect_col %in% names(pleio_data)) {
        stop("Column '", effect_col, "' not found in data")
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    if (show_error_bars && is.null(error_col)) {
        stop("error_col must be specified when show_error_bars = TRUE")
    }

    snp_counts <- pleio_data |>
        group_by(SNPS) |>
        summarise(N_TRAITS = n(), .groups = "drop") |>
        arrange(desc(N_TRAITS)) |>
        slice_head(n = top_n_snps)

    plot_data <- pleio_data |>
        filter(SNPS %in% snp_counts$SNPS) |>
        group_by(SNPS, MAPPED_TRAIT) |>
        summarise(
            effect = max(!!sym(effect_col)),
            .groups = "drop"
        ) |>
        mutate(SNPS = factor(SNPS, levels = rev(snp_counts$SNPS)))

    if (show_error_bars && error_col %in% names(pleio_data)) {
        plot_data <- plot_data |>
            left_join(
                pleio_data |>
                    group_by(SNPS, MAPPED_TRAIT) |>
                    summarise(se = sd(!!sym(error_col), na.rm = TRUE), .groups = "drop"),
                by = c("SNPS", "MAPPED_TRAIT")
            )
    }

    colors <- switch(color_by,
        "snp" = get_pleiotropy_colors("viridis", n = top_n_snps),
        "trait" = get_pleiotropy_colors("okabe_ito", n = length(unique(plot_data$MAPPED_TRAIT))),
        "chromosome" = get_pleiotropy_colors("default", n = 24),
        get_pleiotropy_colors("viridis", n = top_n_snps)
    )

    if (color_by == "snp") {
        p <- ggplot(plot_data, aes(x = SNPS, y = effect, fill = SNPS))
    } else if (color_by == "trait") {
        p <- ggplot(plot_data, aes(x = SNPS, y = effect, fill = MAPPED_TRAIT))
    } else {
        p <- ggplot(plot_data, aes(x = SNPS, y = effect, fill = SNPS))
    }

    p <- p +
        geom_bar(stat = "identity", alpha = 0.8) +
        labs(
            title = "Effect Size Comparison Across Traits",
            subtitle = paste("Top", top_n_snps, "pleiotropic SNPs"),
            x = "SNP",
            y = ifelse(use_log_scale, expression(-log[10](p)), "Effect Size")
        ) +
        theme_pleiotropy_publication() +
        theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)
        )

    if (show_error_bars && "se" %in% names(plot_data)) {
        p <- p + geom_errorbar(aes(ymin = effect - se, ymax = effect + se),
            width = 0.2, alpha = 0.5
        )
    }

    if (color_by == "trait") {
        p <- p + scale_fill_manual(values = colors, name = "Trait")
    } else {
        p <- p + scale_fill_manual(values = colors, guide = "none")
    }

    if (use_log_scale) {
        p <- p + scale_y_continuous(trans = "log10")
    }

    if (flip_coords) {
        p <- p + coord_flip()
    }

    return(p)
}
