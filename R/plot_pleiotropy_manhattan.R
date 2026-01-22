#' Create Manhattan Plot for Pleiotropic SNPs
#'
#' Generates a Manhattan plot highlighting pleiotropic SNPs across chromosomes.
#' The plot shows the distribution of significant associations and can highlight
#' specific SNPs of interest.
#'
#' @param pleio_data A data.frame containing pleiotropy analysis results.
#' @param highlight_snp Character. SNP identifier to highlight (default: NULL).
#' @param title Character. Plot title (default: "Manhattan Plot of Pleiotropic SNPs").
#'
#' @return A ggplot2 object representing the Manhattan plot.
#'
#' @import ggplot2
#' @importFrom dplyr ungroup mutate arrange filter group_by setdiff
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' if (nrow(pleio_results) > 0) {
#'     p <- plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")
#'     print(p)
#' }
#'
#' @import dplyr
#' @export
plot_pleiotropy_manhattan <- function(pleio_data,
                                      highlight_snp = NULL,
                                      title = "Manhattan Plot of Pleiotropic SNPs") {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame or data.table")
    }

    required_cols <- c("SNPS", "CHR_ID", "CHR_POS", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(pleio_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    if (nrow(pleio_data) == 0) {
        stop("Input data is empty")
    }

    plot_data <- pleio_data |>
        mutate(
            CHR_ID = as.character(CHR_ID),
            CHR_POS = as.numeric(CHR_POS),
            CHR_NUM = as.numeric(ifelse(CHR_ID == "X", "23",
                ifelse(CHR_ID == "Y", "24", CHR_ID)
            ))
        ) |>
        filter(!is.na(CHR_NUM), !is.na(CHR_POS)) |>
        arrange(CHR_NUM, CHR_POS)

    # if (nrow(plot_data) == 0) {
    #     stop("No valid chromosomal positions found")
    # }

    chr_lengths <- plot_data |>
        group_by(CHR_ID, CHR_NUM) |>
        summarise(MAX_POS = max(CHR_POS), .groups = "drop") |>
        arrange(CHR_NUM) |>
        mutate(
            CUMPOS_START = cumsum(c(0, MAX_POS[-n()])),
            CUMPOS_END = cumsum(MAX_POS),
            CENTER = CUMPOS_START + (MAX_POS / 2)
        )

    plot_data <- plot_data |>
        left_join(chr_lengths %>% select(CHR_ID, CUMPOS_START), by = "CHR_ID") |>
        mutate(PLOT_POS = CUMPOS_START + CHR_POS)

    colors <- rep(c("#1f77b4", "#ff7f0e"), length.out = nrow(chr_lengths))
    names(colors) <- chr_lengths$CHR_ID

    p <- ggplot(plot_data, aes(
        x = PLOT_POS, y = PVALUE_MLOG,
        color = CHR_ID
    )) +
        geom_point(alpha = 0.7, size = 1.5) +
        scale_color_manual(values = colors) +
        scale_x_continuous(
            labels = chr_lengths$CHR_ID,
            breaks = chr_lengths$CENTER,
            expand = c(0.01, 0)
        ) +
        scale_y_continuous(expand = c(0.02, 0)) +
        labs(
            x = "Chromosome",
            y = expression(-log[10](P)),
            title = title
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, size = 14)
        )

    if (!is.null(highlight_snp) && highlight_snp %in% plot_data$SNPS) {
        highlight_data <- plot_data |>
            filter(SNPS == highlight_snp)

        p <- p +
            geom_point(data = highlight_data, color = "red", size = 3, alpha = 0.8) +
            geom_text(
                data = highlight_data,
                aes(label = SNPS),
                vjust = -0.5, color = "red", size = 3
            )
    }

    return(p)
}
