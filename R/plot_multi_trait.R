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
#' @param suggestive_line Numeric. Suggestive significance threshold (default: 1e-5)
#' @param highlight_snp Character. SNP to highlight across all traits (default: NULL)
#' @param point_size Numeric. Point size (default: 1)
#' @param alpha Numeric. Point transparency (default: 0.6)
#' @param ncol Integer. Number of columns in facet grid (default: 2)
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
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
            CHR_POS = suppressWarnings(as.numeric(CHR_POS)),
            CHR_NUM = suppressWarnings(as.numeric(ifelse(CHR_ID == "X", "23",
                ifelse(CHR_ID == "Y", "24", CHR_ID)
            )))
        ) |>
        filter(!is.na(CHR_NUM), !is.na(CHR_POS)) |>
        arrange(MAPPED_TRAIT, CHR_NUM, CHR_POS)

    if (nrow(plot_data) == 0) {
        warning("No data after filtering")
        return(ggplot() +
            theme_void())
    }

    chr_lengths <- plot_data |>
        group_by(CHR_ID, CHR_NUM) |>
        summarise(MAX_POS = max(CHR_POS), .groups = "drop") |>
        arrange(CHR_NUM) |>
        mutate(
            CUMPOS_START = cumsum(c(0, MAX_POS[-dplyr::n()])),
            CUMPOS_END = cumsum(MAX_POS),
            CENTER = CUMPOS_START + (MAX_POS / 2)
        )

    plot_data <- plot_data |>
        left_join(
            chr_lengths %>% select(CHR_ID, CUMPOS_START),
            by = "CHR_ID"
        ) |>
        mutate(PLOT_POS = CUMPOS_START + CHR_POS)

    chr_levels <- as.character(chr_lengths$CHR_ID)
    colors <- rep(c("#1f77b4", "#ff7f0e"), length.out = length(chr_levels))
    names(colors) <- chr_levels

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
        facet_wrap(~MAPPED_TRAIT, ncol = ncol, scales = "free_y") +
        scale_color_manual(values = colors) +
        scale_x_continuous(
            breaks = chr_lengths$CENTER,
            labels = chr_levels,
            expand = c(0.01, 0)
        ) +
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
        } else {
            warning("highlight_snp '", highlight_snp, "' not found in data")
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
#' @param use_log_scale Logical. Use a log10 scale for the effect sizes. Ignored
#'   when `effect_col` is already a -log10(p) column (default: FALSE)
#' @param show_error_bars Logical. Show error bars (default: FALSE)
#' @param error_col Character. Column holding the standard error of `effect_col`
#'   (default: NULL)
#' @param color_by Character. Color by: "snp", "trait", "chromosome" (default: "snp")
#' @param flip_coords Logical. Flip x and y coordinates (default: TRUE for horizontal bars)
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
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
                                        use_log_scale = FALSE,
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

    if (show_error_bars && !error_col %in% names(pleio_data)) {
        stop("Column '", error_col, "' not found in data")
    }

    if (color_by == "chromosome" && !"CHR_ID" %in% names(pleio_data)) {
        stop("color_by = 'chromosome' requires a 'CHR_ID' column")
    }

    is_neglog_p <- identical(effect_col, "PVALUE_MLOG") ||
        grepl("m?log", effect_col, ignore.case = TRUE)

    snp_counts <- pleio_data |>
        group_by(SNPS) |>
        summarise(N_TRAITS = n(), .groups = "drop") |>
        arrange(desc(N_TRAITS)) |>
        slice_head(n = top_n_snps)

    plot_data <- pleio_data |>
        filter(SNPS %in% snp_counts$SNPS) |>
        group_by(SNPS, MAPPED_TRAIT) |>
        slice_max(!!sym(effect_col), n = 1, with_ties = FALSE) |>
        ungroup() |>
        mutate(
            effect = .data[[effect_col]],
            SNPS = factor(SNPS, levels = rev(snp_counts$SNPS))
        )

    if (show_error_bars) {
        plot_data$se <- suppressWarnings(as.numeric(plot_data[[error_col]]))
    }

    if (color_by == "snp") {
        colors <- get_pleiotropy_colors("viridis", n = top_n_snps)
        p <- ggplot(
            plot_data,
            aes(x = SNPS, y = effect, fill = SNPS, group = MAPPED_TRAIT)
        )
    } else if (color_by == "trait") {
        n_traits <- length(unique(plot_data$MAPPED_TRAIT))
        colors <- get_pleiotropy_colors("okabe_ito", n = n_traits)
        p <- ggplot(
            plot_data,
            aes(x = SNPS, y = effect, fill = MAPPED_TRAIT, group = MAPPED_TRAIT)
        )
    } else {
        chroms <- unique(as.character(plot_data$CHR_ID))
        colors <- rep(c("#1f77b4", "#ff7f0e"), length.out = length(chroms))
        names(colors) <- chroms
        plot_data$CHR_ID <- as.character(plot_data$CHR_ID)
        p <- ggplot(
            plot_data,
            aes(
                x = SNPS, y = effect,
                fill = CHR_ID, group = MAPPED_TRAIT
            )
        )
    }

    dodge <- position_dodge(width = 0.9)

    p <- p +
        geom_col(position = dodge, alpha = 0.8) +
        labs(
            title = "Effect Size Comparison Across Traits",
            subtitle = paste("Top", top_n_snps, "pleiotropic SNPs"),
            x = "SNP",
            y = if (is_neglog_p) expression(-log[10](p)) else "Effect Size"
        ) +
        theme_pleiotropy_publication() +
        theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8)
        )

    if (show_error_bars) {
        p <- p + geom_errorbar(
            aes(ymin = effect - se, ymax = effect + se),
            position = dodge, width = 0.2, alpha = 0.5
        )
    }

    if (color_by == "trait") {
        p <- p + scale_fill_manual(values = colors, name = "Trait")
    } else if (color_by == "chromosome") {
        p <- p + scale_fill_manual(values = colors, name = "Chromosome")
    } else {
        p <- p + scale_fill_manual(values = colors, guide = "none")
    }

    if (use_log_scale && !is_neglog_p) {
        p <- p + scale_y_log10()
    }

    if (flip_coords) {
        p <- p + coord_flip()
    }

    return(p)
}
