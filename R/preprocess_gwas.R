#' Preprocess GWAS Data
#'
#' Filters and cleans GWAS summary statistics for downstream pleiotropy analysis.
#' Removes missing values, applies significance thresholds, and standardizes
#' column names.
#'
#' @param gwas_data A data.frame containing GWAS summary statistics.
#' @param pvalue_threshold Numeric. P-value threshold for filtering
#'   (default: 5e-8).
#' @param columns Character vector. Columns to retain in output
#'   (default: key GWAS columns).
#'
#' @return A filtered and cleaned data.table.
#'
#' @importFrom dplyr filter select all_of setdiff
#'
#' @examples
#' data(gwas_subset)
#' gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
#' head(gwas_clean)
#'
#' @export
preprocess_gwas <- function(gwas_data,
                            pvalue_threshold = 5e-8,
                            columns = c(
                                "SNPS", "MAPPED_TRAIT", "PVALUE_MLOG",
                                "CHR_ID", "CHR_POS"
                            )) {
    if (!is.data.frame(gwas_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (nrow(gwas_data) == 0) {
        stop("Input data is empty")
    }

    names(gwas_data) <- gsub(" ", "_", names(gwas_data))

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(gwas_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    available_cols <- intersect(columns, names(gwas_data))
    if (length(available_cols) == 0) {
        stop("None of the specified columns found in data")
    }

    gwas_data <- gwas_data |>
        filter(!is.na(PVALUE_MLOG)) |>
        filter(PVALUE_MLOG >= -log10(pvalue_threshold)) |>
        select(all_of(available_cols))

    if (nrow(gwas_data) == 0) {
        warning("No data remains after filtering. Consider relaxing p-value threshold.")
    }

    return(gwas_data)
}
