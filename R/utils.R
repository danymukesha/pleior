#' Validate GWAS Data Structure
#'
#' Internal function to validate GWAS data structure and column types.
#'
#' @param gwas_data A data.frame to validate.
#' @param required_cols Character vector of required column names.
#'
#' @return Logical. TRUE if validation passes.
#' @keywords internal
validate_gwas_data <- function(gwas_data, required_cols) {
    if (!is.data.frame(gwas_data)) {
        return(FALSE)
    }

    if (nrow(gwas_data) == 0) {
        return(FALSE)
    }

    missing_cols <- setdiff(required_cols, names(gwas_data))
    if (length(missing_cols) > 0) {
        return(FALSE)
    }

    return(TRUE)
}

#' Format P-values for Display
#'
#' Internal function to format p-values for plotting and output.
#'
#' @param pvalue_mlog Numeric vector of -log10 transformed p-values.
#'
#' @return Character vector of formatted p-values.
#' @keywords internal
format_pvalues <- function(pvalue_mlog) {
    pvalues <- 10^(-pvalue_mlog)
    ifelse(is.na(pvalues), NA_character_,
        ifelse(pvalues < 1e-300, "< 1e-300",
            formatC(pvalues, format = "e", digits = 2)
        )
    )
}

#' Split Multi-Trait Fields into One Row per Trait
#'
#' GWAS Catalog `MAPPED_TRAIT` values may contain several traits separated by
#' commas or semicolons (e.g. "body mass index, obesity"). Counting such a field
#' as a single trait understates pleiotropy, so this helper explodes it into one
#' row per individual trait.
#'
#' @param gwas_data A data.frame.
#' @param trait_col Character. Name of the trait column (default: "MAPPED_TRAIT").
#'
#' @return A data.frame with one row per (association, individual trait).
#' @keywords internal
split_mapped_traits <- function(gwas_data, trait_col = "MAPPED_TRAIT") {
    if (!trait_col %in% names(gwas_data)) {
        return(gwas_data)
    }

    pieces <- strsplit(as.character(gwas_data[[trait_col]]), "[,;]")
    n_pieces <- lengths(pieces)

    out <- gwas_data[rep(seq_len(nrow(gwas_data)), n_pieces), , drop = FALSE]
    out[[trait_col]] <- trimws(unlist(pieces, use.names = FALSE))
    out <- out[!is.na(out[[trait_col]]) & out[[trait_col]] != "", , drop = FALSE]
    rownames(out) <- NULL

    out
}

#' Adjust p-values for Multiple Testing
#'
#' Converts the stored -log10(p-values), applies a multiple-testing correction
#' across all supplied tests, and returns the -log10 transformed adjusted
#' p-values in the same column.
#'
#' @param gwas_data A data.frame.
#' @param method Character. One of "none", "bonferroni" or "fdr".
#' @param pvalue_col Character. Name of the -log10(p) column
#'   (default: "PVALUE_MLOG").
#'
#' @return A data.frame with adjusted -log10(p-values).
#' @importFrom stats p.adjust
#' @keywords internal
adjust_gwas_pvalues <- function(gwas_data, method = "none",
                                pvalue_col = "PVALUE_MLOG") {
    method <- match.arg(method, c("none", "bonferroni", "fdr"))

    if (method == "none" || !pvalue_col %in% names(gwas_data)) {
        return(gwas_data)
    }

    pvalues <- 10^(-gwas_data[[pvalue_col]])
    adjusted <- stats::p.adjust(
        pvalues,
        method = if (method == "fdr") "BH" else "bonferroni"
    )
    gwas_data[[pvalue_col]] <- -log10(adjusted)

    gwas_data
}

#' Distance-Based Clumping of Co-located SNPs
#'
#' A lightweight proxy for LD clumping: within each chromosome, SNPs closer than
#' `clump_kb` kilobases are collapsed and only the lead (most significant) locus
#' of each window is retained. All associations belonging to a retained lead
#' locus are kept, so every trait of a pleiotropic SNP survives. Genotype/LD
#' data are not required, but unlike true LD clumping this method is purely
#' position-based.
#'
#' @param gwas_data A data.frame.
#' @param clump_kb Numeric. Clumping window in kilobases (default: 500).
#' @param chr_col Character. Chromosome column (default: "CHR_ID").
#' @param pos_col Character. Base-pair position column (default: "CHR_POS").
#' @param pvalue_col Character. Significance column, higher is more significant
#'   (default: "PVALUE_MLOG").
#'
#' @return A data.frame retaining the lead SNP of each window.
#' @keywords internal
clump_by_position <- function(gwas_data, clump_kb = 500,
                              chr_col = "CHR_ID", pos_col = "CHR_POS",
                              pvalue_col = "PVALUE_MLOG") {
    if (!all(c(chr_col, pos_col) %in% names(gwas_data))) {
        warning("Clumping requires '", chr_col, "' and '", pos_col,
            "' columns; returning unclumped data",
            call. = FALSE
        )
        return(gwas_data)
    }

    d <- gwas_data
    chr <- as.character(d[[chr_col]])
    pos <- suppressWarnings(as.numeric(d[[pos_col]]))
    sig <- d[[pvalue_col]]
    window <- clump_kb * 1000

    keep <- rep(FALSE, nrow(d))
    keep[is.na(pos)] <- TRUE

    for (this_chr in unique(chr[!is.na(pos)])) {
        idx <- which(!is.na(pos) & chr == this_chr)
        idx <- idx[order(-sig[idx])]

        accepted_pos <- numeric(0)
        for (i in idx) {
            if (length(accepted_pos) == 0 ||
                all(abs(accepted_pos - pos[i]) > window)) {
                accepted_pos <- c(accepted_pos, pos[i])
            }
        }

        # Retain every association of the accepted lead loci so that all traits
        # of a pleiotropic lead SNP are preserved.
        keep[idx[pos[idx] %in% accepted_pos]] <- TRUE
    }

    d <- d[keep, , drop = FALSE]
    rownames(d) <- NULL

    d
}
