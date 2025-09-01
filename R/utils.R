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
  ifelse(pvalues < 1e-300, "< 1e-300",
    formatC(pvalues, format = "e", digits = 2)
  )
}
