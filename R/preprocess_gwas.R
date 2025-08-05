#' Preprocess GWAS Data
#'
#' Filters and cleans GWAS summary statistics for downstream analysis.
#'
#' @param gwas_data A data.table containing GWAS summary statistics.
#' @param pvalue_threshold Numeric. P-value threshold for filtering (default: 5e-8).
#' @param columns Character vector. Columns to retain (default: key GWAS columns).
#' @return A filtered data.table.
#' @importFrom dplyr setdiff filter select
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data, pvalue_threshold = 5e-8)
#' }
#' @export
preprocess_gwas <- function(gwas_data, pvalue_threshold = 5e-8,
                            columns = c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG", "CHR_ID", "CHR_POS")) {
    if (!is.data.frame(gwas_data)) {
        stop("Input must be a data.frame")
    }

    names(gwas_data) <- gsub(" ", "_", names(gwas_data))

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(gwas_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    gwas_data <- gwas_data |>
        filter(!is.na(PVALUE_MLOG)) |>
        filter(PVALUE_MLOG >= -log10(pvalue_threshold))

    gwas_data <- gwas_data |> dplyr::select(all_of(columns))

    return(gwas_data)
}
