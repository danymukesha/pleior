# Prepare Real GWAS Data for pleior Package

# This script loads real GWAS Catalog data and creates a curated subset
# suitable for demonstrating pleiotropy analysis

library(data.table)
library(dplyr)

# Load the full GWAS catalog associations
cat("Loading GWAS catalog associations...\n")
gwas_full <- data.table::fread("data-raw/gwas/gwas_associations.tsv",
    sep = "\t",
    header = TRUE,
    quote = "",
    na.strings = c("NA", "")
)

cat("Original dataset dimensions:", nrow(gwas_full), "rows,", ncol(gwas_full), "columns\n")

# Keep only relevant columns
required_cols <- c(
    "SNPS", "CHR_ID", "CHR_POS", "MAPPED_TRAIT",
    "P-VALUE", "PVALUE_MLOG", "P-VALUE (TEXT)",
    "DISEASE/TRAIT", "REPORTED GENE(S)"
)

gwas_subset <- gwas_full %>%
    select(all_of(intersect(names(.), required_cols)))

# Clean and standardize data
gwas_clean <- gwas_subset %>%
    mutate(
        # Standardize trait names
        MAPPED_TRAIT = ifelse(
            is.na(MAPPED_TRAIT) | MAPPED_TRAIT == "",
            `DISEASE/TRAIT`,
            MAPPED_TRAIT
        ),
        # Ensure PVALUE_MLOG exists
        PVALUE_MLOG = ifelse(
            is.na(PVALUE_MLOG) | PVALUE_MLOG == "",
            -log10(as.numeric(`P-VALUE`)),
            as.numeric(PVALUE_MLOG)
        ),
        # Clean SNP names
        SNPS = gsub("^\\s+", "", SNPS),
        SNPS = gsub("\\s+$", "", SNPS)
    ) %>%
    filter(!is.na(SNPS), !is.na(PVALUE_MLOG), !is.na(MAPPED_TRAIT)) %>%
    filter(PVALUE_MLOG >= 0) %>%
    select(SNPS, MAPPED_TRAIT, PVALUE_MLOG, CHR_ID, CHR_POS)

# Focus on significant associations
gwas_significant <- gwas_clean %>%
    filter(PVALUE_MLOG >= 7) # p < 1e-7

cat("Significant associations:", nrow(gwas_significant), "rows\n")

# Identify traits with most associations
trait_counts <- gwas_significant %>%
    group_by(MAPPED_TRAIT) %>%
    summarise(N_SNPS = n(), .groups = "drop") %>%
    arrange(desc(N_SNPS)) %>%
    head(15)

cat("\nTop traits by number of associations:\n")
print(trait_counts, n = 15)

# Select diverse traits for pleiotropy demonstration
selected_traits <- trait_counts$MAPPED_TRAIT[1:5]
cat("\nSelected traits for pleiotropy analysis:\n")
print(selected_traits)

# Create subset for selected traits
gwas_final <- gwas_clean %>%
    filter(MAPPED_TRAIT %in% selected_traits)

cat("\nFinal subset dimensions:", nrow(gwas_final), "rows\n")

# Count pleiotropic SNPs in the subset
snp_counts <- gwas_final %>%
    group_by(SNPS) %>%
    summarise(N_TRAITS = n(), .groups = "drop") %>%
    filter(N_TRAITS > 1)

cat("Pleiotropic SNPs found:", nrow(snp_counts), "\n")

# Show example pleiotropic SNPs
if (nrow(snp_counts) > 0) {
    cat("\nTop 10 pleiotropic SNPs:\n")
    top_pleiotropic <- snp_counts %>%
        arrange(desc(N_TRAITS)) %>%
        head(10) %>%
        left_join(
            gwas_final %>%
                select(SNPS, MAPPED_TRAIT, PVALUE_MLOG) %>%
                group_by(SNPS, MAPPED_TRAIT) %>%
                summarise(MAX_PVALUE = max(PVALUE_MLOG), .groups = "drop"),
            by = "SNPS"
        )

    print(
        top_pleiotropic %>%
            select(SNPS, N_TRAITS, MAPPED_TRAIT, MAX_PVALUE),
        n = 10
    )
}

# Save as .rda file for package
output_path <- "data/gwas_subset.rda"
save(gwas_final, file = output_path)
cat("\nSaved to:", output_path, "\n")

# Also save as .rds for flexibility
output_rds <- "data/gwas_subset.rds"
saveRDS(gwas_final, file = output_rds)
cat("Also saved to:", output_rds, "\n")

usethis::use_data(gwas_subset, overwrite = TRUE)

# Summary statistics
cat("\n=== DATA PREPARATION SUMMARY ===\n")
cat("Original rows:", nrow(gwas_full), "\n")
cat("After cleaning:", nrow(gwas_clean), "\n")
cat("Final subset:", nrow(gwas_final), "\n")
cat("Unique SNPs:", length(unique(gwas_final$SNPS)), "\n")
cat("Unique traits:", length(unique(gwas_final$MAPPED_TRAIT)), "\n")
cat("Pleiotropic SNPs:", nrow(snp_counts), "\n")
cat("\nTraits included:\n")
print(unique(gwas_final$MAPPED_TRAIT))
