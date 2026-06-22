#' pleior: Pleiotropy Analysis for GWAS Data
#'
#' The pleior package provides tools for identifying and analyzing pleiotropic genetic variants
#' in genome-wide association studies (GWAS). It includes functions for data loading,
#' preprocessing, pleiotropy detection, and visualization, enabling researchers to explore
#' shared genetic mechanisms across complex traits.
#'
#' @name pleior-package
#' @author Dany Mukesha
#' @references
#' Watanabe, K., et al. (2019). A global overview of pleiotropy and genetic architecture
#' in complex traits. *Nature Genetics*, 51, 1339–1348.
#' @keywords GWAS pleiotropy genetics
#' @importFrom grDevices colorRampPalette
#' @importFrom stats density dist hclust median quantile sd
"_PACKAGE"

utils::globalVariables(c(
    "CHR_ID", "CHR_POS", "CHR_NUM", "CUMPOS_START", "MAX_POS",
    "PLOT_POS", "PVALUE_MLOG", "SNPS", "MAPPED_TRAIT", "N_TRAITS",
    "N_SNPS", "TRAIT", "CENTER", "value", "total", "trait1", "trait2",
    "shared_snps", "weight", "n_snps", "type", "n_traits", "name",
    "MAX_PVALUE", "distance_kb", "is_target", "CHR", "START", "END",
    "mid_pos", "effect", "se", "x", "y", "r", "color", "density"
))
