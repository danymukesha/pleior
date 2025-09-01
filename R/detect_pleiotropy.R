#' Detect Pleiotropic SNPs
#'
#' Identifies SNPs associated with multiple traits in GWAS data. A SNP is
#' considered pleiotropic if it shows significant associations with two or
#' more distinct traits.
#'
#' @param gwas_data A data.frame containing preprocessed GWAS summary statistics.
#' @param traits Character vector. Specific traits to analyze for pleiotropy
#'   (default: NULL uses all traits).
#' @param pvalue_threshold Numeric. P-value threshold for significance
#'   (default: 5e-8).
#'
#' @return A data.table with pleiotropic SNPs, their associated traits,
#'   and significance levels.
#'
#' @importFrom dplyr group_by summarise filter left_join n setdiff
#' @importFrom stringr str_detect
#'
#' @examples
#' data(gwas_subset)
#' pleio_results <- detect_pleiotropy(gwas_subset)
#' head(pleio_results)
#'
#' # Analyze specific traits
#' specific_traits <- c("Alzheimer disease", "myocardial infarction")
#' pleio_specific <- detect_pleiotropy(gwas_subset, traits = specific_traits)
#'
#' @export
detect_pleiotropy <- function(gwas_data, traits = NULL, pvalue_threshold = 5e-8) {
    if (!is.data.frame(gwas_data)) {
        stop("Input must be a data.frame or data.table")
    }

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(gwas_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    if (nrow(gwas_data) == 0) {
        stop("Input data is empty")
    }

    if (!is.null(traits)) {
        trait_pattern <- paste(traits, collapse = "|")
        gwas_data <- gwas_data |>
            filter(str_detect(MAPPED_TRAIT, trait_pattern))

        if (nrow(gwas_data) == 0) {
            stop("No data found for specified traits")
        }
    }

    pleio_table <- gwas_data |>
        group_by(SNPS, MAPPED_TRAIT) |>
        summarise(
            TRAIT = MAPPED_TRAIT[which.max(PVALUE_MLOG)],
            MAX_PVALUE_MLOG = max(PVALUE_MLOG),
            .groups = "drop"
        ) |>
        group_by(SNPS) |>
        summarise(
            N_TRAITS = n(),
            TRAITS = paste(unique(TRAIT), collapse = ";"),
            .groups = "drop"
        ) |>
        filter(N_TRAITS > 1)

    if (nrow(pleio_table) == 0) {
        warning("No pleiotropic SNPs found with current parameters")
        return(data.table::data.table())
    }

    pleio_results <- pleio_table |>
        left_join(gwas_data, by = "SNPS") |>
        filter(PVALUE_MLOG >= -log10(pvalue_threshold))

    return(data.table::as.data.table(pleio_results))
}
