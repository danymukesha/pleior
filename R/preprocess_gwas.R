#' Preprocess GWAS Data
#'
#' Filters and cleans GWAS summary statistics for downstream analysis.
#'
#' @param gwas_data A data.table containing GWAS summary statistics.
#' @param pvalue_threshold Numeric. P-value threshold for filtering (default: 5e-8).
#' @param columns Character vector. Columns to retain (default: key GWAS columns).
#' @return A filtered data.table.
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data, pvalue_threshold = 5e-8)
#' }
#' @export
preprocess_gwas <- function(gwas_data, pvalue_threshold = 5e-8,
                            columns = c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG", "CHR_ID", "CHR_POS")) {
    if (!inherits(gwas_data, "data.table")) {
        stop("Input must be a data.table")
    }
    print("he")
    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    if (!all(required_cols %in% names(gwas_data))) {
        stop("Required columns missing: ", paste(required_cols[!required_cols %in% names(gwas_data)], collapse = ", "))
    }
    gwas_data <- gwas_data[!is.na(PVALUE_MLOG) & PVALUE_MLOG >= -log10(pvalue_threshold), ]
    gwas_data[, .SD, .SDcols = columns]
}
