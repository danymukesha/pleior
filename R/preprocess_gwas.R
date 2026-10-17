#' Preprocess GWAS Data
#'
#' Filters and cleans GWAS summary statistics for downstream pleiotropy analysis.
#' Removes missing values, applies significance thresholds, optionally adjusts
#' for multiple testing, and standardizes column names.
#'
#' The -log10(p-value) column (`PVALUE_MLOG`) is derived automatically from a
#' raw p-value column (`P-VALUE`, `PVALUE` or `P`) when it is not present.
#'
#' @param gwas_data A data.frame containing GWAS summary statistics.
#' @param pvalue_threshold Numeric. P-value threshold for filtering
#'   (default: 5e-8).
#' @param columns Character vector. Columns to retain in output
#'   (default: key GWAS columns).
#' @param adjust_method Character. Multiple-testing correction applied before
#'   thresholding: "none" (default), "bonferroni" or "fdr".
#'
#' @return A filtered and cleaned data.frame (of the same class as the input).
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
                            ),
                            adjust_method = c("none", "bonferroni", "fdr")) {
    if (!is.data.frame(gwas_data)) {
        stop("Input must be a data.frame or data.table")
    }

    if (nrow(gwas_data) == 0) {
        stop("Input data is empty")
    }

    adjust_method <- match.arg(adjust_method)

    data <- as.data.frame(gwas_data)
    names(data) <- gsub(" ", "_", names(data))
    columns <- gsub(" ", "_", columns)

    if (!"PVALUE_MLOG" %in% names(data)) {
        p_candidates <- intersect(c("P-VALUE", "P_VALUE", "PVALUE", "P"), names(data))
        if (length(p_candidates) > 0) {
            raw_p <- suppressWarnings(as.numeric(data[[p_candidates[1]]]))
            data$PVALUE_MLOG <- -log10(pmax(raw_p, 1e-300, na.rm = FALSE))
        }
    }

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    available_cols <- intersect(columns, names(data))
    if (length(available_cols) == 0) {
        stop("None of the specified columns found in data")
    }

    data$PVALUE_MLOG <- suppressWarnings(as.numeric(data$PVALUE_MLOG))
    data <- data[!is.na(data$PVALUE_MLOG), , drop = FALSE]

    data <- adjust_gwas_pvalues(data, method = adjust_method)

    data <- data[data$PVALUE_MLOG >= -log10(pvalue_threshold), , drop = FALSE]
    data <- data[, available_cols, drop = FALSE]
    rownames(data) <- NULL

    if (nrow(data) == 0) {
        warning("No data remains after filtering. Consider relaxing p-value threshold.")
    }

    if (inherits(gwas_data, "data.table")) {
        data <- data.table::as.data.table(data)
    }

    data
}
