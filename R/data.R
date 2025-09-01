#' Example GWAS subset data
#'
#' A subset of GWAS summary statistics containing SNPs and their associations
#' with multiple traits, used for demonstrating pleiotropy analysis.
#'
#' @format A data.table with 10 rows and 5 columns:
#' \describe{
#'   \item{SNPS}{Character. SNP identifiers (rs numbers)}
#'   \item{MAPPED_TRAIT}{Character. Associated trait or phenotype}
#'   \item{PVALUE_MLOG}{Numeric. -log10 transformed p-values}
#'   \item{CHR_ID}{Character. Chromosome identifier}
#'   \item{CHR_POS}{Character. Chromosomal position}
#' }
#'
#' @source Simulated data based on real GWAS catalog associations
#' @examples
#' data(gwas_subset)
#' head(gwas_subset)
"gwas_subset"
