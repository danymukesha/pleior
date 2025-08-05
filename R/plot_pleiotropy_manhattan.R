#' Plot Manhattan Plot for Pleiotropic SNPs
#'
#' Creates a Manhattan plot highlighting pleiotropic SNPs across chromosomes.
#'
#' @param pleio_data A data.table containing pleiotropy analysis results.
#' @param highlight_snp Character. Optional SNP to highlight (default: NULL).
#' @return A ggplot2 object.
#' @import ggplot2
#' @importFrom dplyr ungroup mutate arrange filter group_by
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data)
#' pleio_results <- detect_pleiotropy(gwas_clean)
#' plot_pleiotropy_manhattan(pleio_results, highlight_snp = "rs814573")
#' }
#' @export
plot_pleiotropy_manhattan <- function(pleio_data, highlight_snp = NULL) {
    if (!is.data.frame(pleio_data)) {
        stop("Input must be a data.frame")
    }

    # Check if required columns are present
    required_cols <- c("SNPS", "CHR_ID", "CHR_POS", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(pleio_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    # Prepare the plot data
    plot_data <- pleio_data |>
        mutate(CHR_POS = as.numeric(CHR_POS)) |>
        group_by(CHR_ID) |>
        arrange(CHR_POS) |>
        mutate(CHR_CUMPOS = CHR_POS + cumsum(c(0, diff(CHR_POS, lag = 1, differences = 1))[1:n()])) |>
        ungroup()

    # Calculate x-axis positioning for chromosomes
    axis_df <- plot_data |>
        group_by(CHR_ID) |>
        summarise(CENTER = mean(CHR_CUMPOS), .groups = "drop")

    # Create the basic Manhattan plot
    p <- ggplot(plot_data, aes(x = CHR_CUMPOS, y = PVALUE_MLOG, color = as.factor(CHR_ID))) +
        geom_point() +
        scale_x_continuous(label = axis_df$CHR_ID, breaks = axis_df$CENTER) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Chromosome", y = "-log10(P-value)", title = "Manhattan Plot of Pleiotropic SNPs") +
        theme_minimal() +
        theme(legend.position = "none")

    # Highlight the selected SNP if provided
    if (!is.null(highlight_snp)) {
        highlight_data <- plot_data |>
            filter(SNPS == highlight_snp)

        if (nrow(highlight_data) > 0) {
            p <- p + geom_point(data = highlight_data, color = "red", size = 3)
        }
    }

    return(p)
}
