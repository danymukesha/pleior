#' Detect Pleiotropic SNPs
#'
#' Identifies SNPs associated with multiple traits in GWAS data. A SNP is
#' considered pleiotropic if it shows significant associations with two or
#' more distinct traits after applying the significance threshold.
#'
#' The significance filter is applied *before* traits are counted, so every
#' reported pleiotropic SNP has at least two associations passing the
#' threshold. Multi-trait `MAPPED_TRAIT` fields (comma- or semicolon-separated)
#' are expanded into individual traits.
#'
#' @param gwas_data A data.frame containing GWAS summary statistics.
#' @param traits Character vector. Specific traits to analyze for pleiotropy
#'   (default: NULL uses all traits). Matching is exact against individual
#'   traits after splitting multi-trait fields.
#' @param pvalue_threshold Numeric. P-value threshold for significance
#'   (default: 5e-8).
#' @param adjust_method Character. Multiple-testing correction applied to all
#'   associations before thresholding: "none" (default), "bonferroni" or "fdr".
#' @param clump_kb Numeric or NULL. If supplied, co-located SNPs within this
#'   distance (kilobases) are clumped to their most significant representative.
#'   This is a position-based proxy for LD clumping (default: NULL, no clumping).
#'
#' @return A data.frame with pleiotropic SNPs, one row per SNP-trait
#'   association, including `N_TRAITS` and `TRAITS` summary columns.
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
detect_pleiotropy <- function(gwas_data, traits = NULL, pvalue_threshold = 5e-8,
                              adjust_method = c("none", "bonferroni", "fdr"),
                              clump_kb = NULL) {
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

    adjust_method <- match.arg(adjust_method)

    data <- as.data.frame(gwas_data)
    data$PVALUE_MLOG <- suppressWarnings(as.numeric(data$PVALUE_MLOG))
    data <- data[!is.na(data$PVALUE_MLOG), , drop = FALSE]

    data <- split_mapped_traits(data)

    data <- adjust_gwas_pvalues(data, method = adjust_method)

    data <- data[data$PVALUE_MLOG >= -log10(pvalue_threshold), , drop = FALSE]

    if (!is.null(traits)) {
        data <- data[data$MAPPED_TRAIT %in% traits, , drop = FALSE]

        if (nrow(data) == 0) {
            stop("No data found for specified traits")
        }
    }

    if (nrow(data) == 0) {
        warning("No pleiotropic SNPs found with current parameters")
        return(dplyr::tibble())
    }

    if (!is.null(clump_kb)) {
        data <- clump_by_position(data, clump_kb = clump_kb)
    }

    snp_trait <- data |>
        group_by(SNPS, MAPPED_TRAIT) |>
        slice_max(PVALUE_MLOG, n = 1, with_ties = FALSE) |>
        ungroup()

    trait_summary <- snp_trait |>
        group_by(SNPS) |>
        summarise(
            N_TRAITS = n(),
            TRAITS = paste(sort(unique(MAPPED_TRAIT)), collapse = ";"),
            .groups = "drop"
        ) |>
        filter(N_TRAITS > 1)

    if (nrow(trait_summary) == 0) {
        warning("No pleiotropic SNPs found with current parameters")
        return(dplyr::tibble())
    }

    result <- snp_trait |>
        filter(SNPS %in% trait_summary$SNPS) |>
        left_join(trait_summary, by = "SNPS") |>
        arrange(desc(N_TRAITS), desc(PVALUE_MLOG), SNPS)

    as.data.frame(result)
}
