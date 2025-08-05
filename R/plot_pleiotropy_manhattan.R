#' Plot Manhattan Plot for Pleiotropic SNPs
#'
#' Creates a Manhattan plot highlighting pleiotropic SNPs across chromosomes.
#'
#' @param pleio_data A data.table containing pleiotropy analysis results.
#' @param highlight_snp Character. Optional SNP to highlight (default: NULL).
#' @return A ggplot2 object.
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data)
#' pleio_results <- detect_pleiotropy(gwas_clean)
#' plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")
#' }
#' @export
plot_pleiotropy_manhattan <- function(pleio_data, highlight_snp = NULL) {
    if (!inherits(pleio_data, "data.table")) {
        stop("Input must be a data.table")
    }
    required_cols <- c("SNPS", "CHR_ID", "CHR_POS", "PVALUE_MLOG")
    if (!all(required_cols %in% names(pleio_data))) {
        stop("Required columns missing: ", paste(required_cols[!required_cols %in% names(pleio_data)], collapse = ", "))
    }

    plot_data <- pleio_data[, .(CHR_ID, CHR_POS = as.numeric(CHR_POS), PVALUE_MLOG, SNPS)]
    plot_data[, CHR_CUMPOS := CHR_POS + cumsum(as.numeric(c(0, max(CHR_POS, na.rm = TRUE)))[1:.N]), by = CHR_ID]

    axis_df <- plot_data[, .(CENTER = mean(CHR_CUMPOS)), by = CHR_ID]

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = CHR_CUMPOS, y = PVALUE_MLOG, color = as.factor(CHR_ID))) +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(label = axis_df$CHR_ID, breaks = axis_df$CENTER) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(x = "Chromosome", y = "-log10(P-value)", title = "Manhattan Plot of Pleiotropic SNPs") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")

    if (!is.null(highlight_snp)) {
        highlight_data <- plot_data[SNPS == highlight_snp]
        if (nrow(highlight_data) > 0) {
            p <- p + ggplot2::geom_point(data = highlight_data, color = "red", size = 3)
        }
    }

    return(p)
}
