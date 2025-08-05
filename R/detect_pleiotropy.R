#' Detect Pleiotropic SNPs
#'
#' Identifies SNPs associated with multiple traits in GWAS data.
#'
#' @param gwas_data A data.table containing preprocessed GWAS summary statistics.
#' @param traits Character vector. Traits to analyze for pleiotropy (default: NULL, uses all traits).
#' @param pvalue_threshold Numeric. P-value threshold for significance (default: 5e-8).
#' @return A data.table with pleiotropic SNPs, their associated traits, and significance levels.
#' @importFrom dplyr group_by summarise filter left_join n
#' @importFrom stringr str_detect
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data)
#' pleio_results <- detect_pleiotropy(gwas_clean, traits = c("Alzheimer disease", "myocardial infarction"))
#' }
#' @export
detect_pleiotropy <- function(gwas_data, traits = NULL, pvalue_threshold = 5e-8) {
    if (!is.data.frame(gwas_data)) {
        stop("Input must be a data.frame")
    }

    required_cols <- c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG")
    missing_cols <- setdiff(required_cols, names(gwas_data))
    if (length(missing_cols) > 0) {
        stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
    }

    if (!is.null(traits)) {
        gwas_data <- gwas_data |>
            filter(str_detect(MAPPED_TRAIT, paste(traits, collapse = "|")))
    }

    # get the max PVALUE_MLOG for each SNP
    pleio_table <- gwas_data |>
        group_by(SNPS, MAPPED_TRAIT) |>
        summarise(
            TRAIT = MAPPED_TRAIT[which.max(PVALUE_MLOG)],
            PVALUE_MLOG = max(PVALUE_MLOG), .groups = "drop"
        )

    # count the number of traits associated with each SNP and collapse them
    pleio_table <- pleio_table |>
        group_by(SNPS) |>
        summarise(
            N_TRAITS = n(),
            TRAITS = paste(TRAIT, collapse = ";"), .groups = "drop"
        ) |>
        filter(N_TRAITS > 1) # i only keep SNPs associated with multiple traits

    # merging back with the original data for further details
    pleio_results <- pleio_table |>
        left_join(gwas_data, by = "SNPS") |>
        filter(PVALUE_MLOG >= -log10(pvalue_threshold))

    return(pleio_results)
}
