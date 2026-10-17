#' GWAS Summary Statistics Subset
#'
#' A curated subset of GWAS Catalog summary statistics prepared for pleiotropy
#' analysis.
#'
#' This dataset contains genome-wide association study (GWAS) summary
#' statistics from the GWAS Catalog (https://www.ebi.ac.uk/gwas/), reduced to
#' a small, reproducible set that demonstrates pleiotropy detection and
#' visualization. It is the example dataset shipped with `pleior`.
#'
#' @format A data.table with 2498 rows and 5 columns:
#' \describe{
#'   \item{SNPS}{Character. SNP identifiers (rs numbers)}
#'   \item{MAPPED_TRAIT}{Character. Standardized trait or phenotype name}
#'   \item{PVALUE_MLOG}{Numeric. -log10 transformed p-values for significance}
#'   \item{CHR_ID}{Character. Chromosome identifier (1-22, X, Y)}
#'   \item{CHR_POS}{Numeric. Chromosomal position in base pairs}
#' }
#'
#' @source
#' GWAS Catalog associations file (2024 release), downloaded from:
#' \url{https://www.ebi.ac.uk/gwas/api/search/downloads/alternative}
#'
#' The subset was created by \code{data-raw/create_gwas_subset.R}:
#' 1. Loading the complete GWAS Catalog associations file
#' 2. Selecting 15 traits relevant to complex disease and metabolism
#' 3. Collapsing duplicate study entries to one row per SNP-trait pair
#'    (keeping the most significant association)
#' 4. Retaining all associations of the strongest pleiotropic SNPs, then
#'    filling with a trait-balanced random sample to reach ~2500 rows
#' 5. Ensuring representation across multiple chromosomes
#'
#' @details
#' The dataset contains 2498 associations covering 1835 distinct SNPs across
#' 15 traits. Of these SNPs, 128 are pleiotropic (associated with two or more
#' of the 15 traits); the strongest, \code{rs429358} (in the APOE region),
#' maps to 12 traits. \code{rs814573}, used as a running example in the
#' README, is associated with myocardial infarction and C-reactive protein
#' measurement in this subset.
#'
#' The 15 traits are:
#' \itemize{
#'   \item body mass index
#'   \item body height
#'   \item Alzheimer disease
#'   \item myocardial infarction
#'   \item coronary artery disease
#'   \item type 2 diabetes mellitus
#'   \item systolic blood pressure
#'   \item diastolic blood pressure
#'   \item low density lipoprotein cholesterol measurement
#'   \item high density lipoprotein cholesterol measurement
#'   \item triglyceride measurement
#'   \item total cholesterol measurement
#'   \item C-reactive protein measurement
#'   \item asthma
#'   \item schizophrenia
#' }
#'
#' The data is suitable for:
#' \itemize{
#'   \item Demonstrating pleiotropy detection (SNPs associated with multiple traits)
#'   \item Testing visualization functions
#'   \item Teaching GWAS analysis workflows
#'   \item Benchmarking analysis methods
#' }
#'
#' @references
#' GWAS Catalog Team. (2024). GWAS Catalog: a knowledgebase for genome-wide
#' association studies. *Nucleic Acids Research*, 51(D1), D944-D949.
#' \url{https://doi.org/10.1093/nar/gkab1009}
#'
#' @examples
#' # Load the data
#' data(gwas_subset)
#' head(gwas_subset)
#'
#' # Basic statistics
#' n_unique_snps <- length(unique(gwas_subset$SNPS))
#' n_unique_traits <- length(unique(gwas_subset$MAPPED_TRAIT))
#' cat("Unique SNPs:", n_unique_snps, "\n")
#' cat("Unique traits:", n_unique_traits, "\n")
#'
#' # Preprocess and detect pleiotropy
#' gwas_clean <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
#' pleio_results <- detect_pleiotropy(gwas_clean)
#' head(pleio_results)
"gwas_subset"
