#' Load GWAS Summary Statistics
#'
#' Loads GWAS summary statistics from a file or URL, supporting various formats.
#' @param file_path Character. Path to the GWAS summary statistics file (e.g., TSV, CSV).
#' @param sep Character. Separator used in the file (default: "\t").
#' @param header Logical. Whether the file has a header row (default: TRUE).
#' @param quote Character. Quote character used in the file (default: "").
#' @return A data.table containing the GWAS data.
#' @examples #'
#' \dontrun{
#' gwas_data <- load_gwas_data("gwas_associations.tsv")
#' }
#' @export
load_gwas_data <- function(file_path, sep = "\t",
                           header = TRUE,
                           quote = "") {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }
    data.table::fread(file_path, sep = sep, header = header, quote = quote)
}
