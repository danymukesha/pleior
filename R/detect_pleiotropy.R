#' Detect Pleiotropic SNPs
#'
#' Identifies SNPs associated with multiple traits in GWAS data.
#'
#' @param gwas_data A data.table containing preprocessed GWAS summary statistics.
#' @param traits Character vector. Traits to analyze for pleiotropy (default: NULL, uses all traits).
#' @param pvalue_threshold Numeric. P-value threshold for significance (default: 5e-8).
#' @return A data.table with pleiotropic SNPs, their associated traits, and significance levels.
#' @examples
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' gwas_clean <- preprocess_gwas(gwas_data)
#' pleio_results <- detect_pleiotropy(gwas_clean, traits = c("Alzheimer disease", "myocardial infarction"))
#' }
#' @export
detect_pleiotropy <- function(gwas_data, traits = NULL, pvalue_threshold = 5e-8) {
    if (!inherits(gwas_data, "data.table")) {
        stop("Input must be a data.table")
    }
    if (!all(c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG") %in% names(gwas_data))) {
        stop("Required columns missing: SNPS, MAPPED_TRAIT, PVALUE_MLOG")
    }

    if (!is.null(traits)) {
        gwas_data <- gwas_data[stringr::str_detect(MAPPED_TRAIT, paste(traits, collapse = "|")), ]
    }

    pleio_table <- gwas_data[, .(TRAIT = MAPPED_TRAIT, PVALUE_MLOG = max(PVALUE_MLOG)), by = SNPS]
    pleio_table <- pleio_table[, .(N_TRAITS = .N, TRAITS = paste(TRAIT, collapse = ";")), by = SNPS]
    pleio_table <- pleio_table[N_TRAITS > 1, ]

    pleio_results <- merge(pleio_table, gwas_data, by = "SNPS", allow.cartesian = TRUE)
    pleio_results[PVALUE_MLOG >= -log10(pvalue_threshold), ]
}
