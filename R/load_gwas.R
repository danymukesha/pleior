#' Load GWAS Summary Statistics
#'
#' Loads GWAS summary statistics from a file, supporting various formats.
#' The function automatically detects common separators and handles
#' standard GWAS file formats.
#'
#' @param file_path Character. Path to the GWAS summary statistics file.
#' @param sep Character. Field separator (default: "\\t").
#' @param header Logical. Whether file contains header row (default: TRUE).
#' @param quote Character. Quote character (default: "").
#'
#' @return A data.table containing the GWAS summary statistics.
#'
#' @importFrom data.table fread
#'
#' @examples
#' \donttest{
#' # Load example data
#' file_path <- system.file("extdata", "example_gwas.tsv", package = "pleior")
#' if (file.exists(file_path)) {
#'     gwas_data <- load_gwas_data(file_path)
#'     head(gwas_data)
#' }
#' }
#'
#' @export
load_gwas_data <- function(file_path, sep = "\t", header = TRUE, quote = "") {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    tryCatch(
        {
            data.table::fread(file_path, sep = sep, header = header, quote = quote)
        },
        error = function(e) {
            stop("Error reading file: ", e$message)
        }
    )
}
