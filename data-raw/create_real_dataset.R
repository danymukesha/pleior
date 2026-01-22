#!/usr/bin/env Rscript

# Create Real GWAS Subset for pleior Package
# This script processes GWAS Catalog data and creates a curated dataset

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

cat("=== Creating Real GWAS Dataset for pleior ===\n\n")

# Define file paths
input_file <- "data-raw/gwas/gwas_associations.tsv"
output_file <- "data/gwas_subset.rda"
output_rds <- "data/gwas_subset.rds"

# Check if input exists
if (!file.exists(input_file)) {
    stop(
        "Input file not found: ", input_file,
        "\nPlease download GWAS catalog data first."
    )
}

cat("Step 1: Loading GWAS catalog associations...\n")

# Load the GWAS catalog file
# Using fread for efficient loading of large files
gwas_full <- fread(
    input_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    na.strings = c("NA", "", " ", "."),
    select = c(
        "SNPS", "CHR_ID", "CHR_POS", "MAPPED_TRAIT",
        "P-VALUE", "PVALUE_MLOG", "DISEASE/TRAIT",
        "REPORTED GENE(S)"
    )
)

cat("Loaded:", nrow(gwas_full), "rows,", ncol(gwas_full), "columns\n")

# Clean and standardize the data
cat("Step 2: Cleaning and standardizing data...\n")

gwas_clean <- gwas_full %>%
    mutate(
        # Use MAPPED_TRAIT if available, otherwise use DISEASE/TRAIT
        MAPPED_TRAIT = ifelse(
            is.na(MAPPED_TRAIT) | MAPPED_TRAIT == "",
            `DISEASE/TRAIT`,
            MAPPED_TRAIT
        ),
        # Calculate PVALUE_MLOG if missing
        PVALUE_MLOG = ifelse(
            is.na(PVALUE_MLOG) | PVALUE_MLOG == "",
            -log10(pmax(as.numeric(`P-VALUE`), 1e-300)),
            as.numeric(PVALUE_MLOG)
        ),
        # Clean SNP identifiers
        SNPS = gsub("\\s+$", "", SNPS),
        SNPS = gsub("^\\s+", "", SNPS)
    ) %>%
    filter(
        !is.na(SNPS),
        !is.na(PVALUE_MLOG),
        !is.na(MAPPED_TRAIT),
        PVALUE_MLOG >= 0,
        PVALUE_MLOG < Inf
    ) %>%
    select(SNPS, MAPPED_TRAIT, PVALUE_MLOG, CHR_ID, CHR_POS)

cat("After cleaning:", nrow(gwas_clean), "rows\n")

# Filter for significant associations (p < 1e-5 to reduce size)
cat("Step 3: Filtering for significant associations (p < 1e-5)...\n")

gwas_significant <- gwas_clean %>%
    filter(PVALUE_MLOG >= 5) # -log10(1e-5) = 5

cat("Significant associations:", nrow(gwas_significant), "rows\n")

# Identify top traits by number of associations
cat("Step 4: Selecting diverse traits...\n")

trait_counts <- gwas_significant %>%
    group_by(MAPPED_TRAIT) %>%
    summarise(N_SNPS = n(), .groups = "drop") %>%
    arrange(desc(N_SNPS))

cat("Found", nrow(trait_counts), "unique traits\n")

# Select top 15 traits for diversity and coverage
top_traits <- trait_counts$MAPPED_TRAIT[1:15]

cat("Selected top 15 traits:\n")
for (i in 1:length(top_traits)) {
    cat(sprintf(
        "  %d. %s (%d SNPs)\n", i, top_traits[i],
        trait_counts$N_SNPS[i]
    ))
}

# Filter for selected traits
gwas_filtered <- gwas_significant %>%
    filter(MAPPED_TRAIT %in% top_traits)

cat("\nAfter trait filtering:", nrow(gwas_filtered), "rows\n")

# Sample to create manageable dataset (aim for ~2000-3000 rows)
target_size <- 2500
if (nrow(gwas_filtered) > target_size) {
    cat("Step 5: Sampling to create manageable dataset (target:", target_size, "rows)...\n")

    # Sample proportionally from each trait
    gwas_final <- gwas_filtered %>%
        group_by(MAPPED_TRAIT) %>%
        sample_frac(min(1, target_size / n())) %>%
        ungroup()
} else {
    gwas_final <- gwas_filtered
}

cat("Final dataset size:", nrow(gwas_final), "rows\n")

# Analyze pleiotropy in the final dataset
cat("\nStep 6: Analyzing pleiotropy...\n")

snp_counts <- gwas_final %>%
    group_by(SNPS) %>%
    summarise(N_TRAITS = n(), .groups = "drop") %>%
    filter(N_TRAITS > 1) %>%
    arrange(desc(N_TRAITS))

cat("Pleiotropic SNPs found:", nrow(snp_counts), "\n")

if (nrow(snp_counts) > 0) {
    cat("\nTop 10 pleiotropic SNPs:\n")
    for (i in 1:min(10, nrow(snp_counts))) {
        snp_id <- snp_counts$SNPS[i]
        n_traits <- snp_counts$N_TRAITS[i]

        # Get traits for this SNP
        traits_list <- gwas_final %>%
            filter(SNPS == snp_id) %>%
            pull(MAPPED_TRAIT) %>%
            unique() %>%
            paste(collapse = ", ")

        cat(sprintf("  %s: %d traits - %s\n", snp_id, n_traits, traits_list))
    }
}

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Original rows:", nrow(gwas_full), "\n")
cat("After cleaning:", nrow(gwas_clean), "\n")
cat("Significant (p<1e-5):", nrow(gwas_significant), "\n")
cat("After trait filtering:", nrow(gwas_filtered), "\n")
cat("Final dataset:", nrow(gwas_final), "\n")
cat("Unique SNPs:", length(unique(gwas_final$SNPS)), "\n")
cat("Unique traits:", length(unique(gwas_final$MAPPED_TRAIT)), "\n")
cat("Pleiotropic SNPs:", nrow(snp_counts), "\n")
cat("\nTraits in final dataset:\n")
print(unique(gwas_final$MAPPED_TRAIT))

# Save the data
cat("\nStep 7: Saving data files...\n")

# Save as .rda for package data() function
save(gwas_final, file = output_file)
cat("Saved:", output_file, "\n")

# Also save as .rds for flexibility
saveRDS(gwas_final, file = output_rds)
cat("Also saved:", output_rds, "\n")

usethis::use_data(gwas_final, overwrite = TRUE)

cat("\n=== SUCCESS ===\n")
cat("Real GWAS dataset created successfully!\n")
