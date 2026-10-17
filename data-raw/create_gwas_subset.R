#!/usr/bin/env Rscript

# Curate a small, honest pleiotropy demonstration dataset from the full
# GWAS Catalog dump shipped in data/gwas_subset.rda.
#
# Usage:
#   Rscript data-raw/create_gwas_subset.R
#
# The script expects the full GWAS Catalog associations (as distributed in the
# previous data/gwas_subset.rda) and writes:
#   - data/gwas_subset.rda         (object: gwas_subset)
#   - inst/extdata/gwas_subset.rds
#   - inst/extdata/example_gwas.tsv

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
input_rda <- if (length(args) >= 1) args[1] else "data/gwas_subset.rda"

load(input_rda) # provides `gwas_subset`
d <- as.data.frame(gwas_subset)

traits15 <- c(
    "body mass index",
    "body height",
    "Alzheimer disease",
    "myocardial infarction",
    "coronary artery disease",
    "type 2 diabetes mellitus",
    "systolic blood pressure",
    "diastolic blood pressure",
    "low density lipoprotein cholesterol measurement",
    "high density lipoprotein cholesterol measurement",
    "triglyceride measurement",
    "total cholesterol measurement",
    "C-reactive protein measurement",
    "asthma",
    "schizophrenia"
)

sub <- d[
    d$MAPPED_TRAIT %in% traits15 &
        !is.na(d$SNPS) & !is.na(d$MAPPED_TRAIT) &
        !is.na(d$PVALUE_MLOG) & !is.na(d$CHR_ID) & !is.na(d$CHR_POS),
    c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG", "CHR_ID", "CHR_POS")
]

sub$CHR_ID <- as.character(sub$CHR_ID)
sub$CHR_POS <- suppressWarnings(as.numeric(sub$CHR_POS))
sub$PVALUE_MLOG <- suppressWarnings(as.numeric(sub$PVALUE_MLOG))

# Collapse duplicate study entries: one row per (SNP, trait), most significant.
sub <- sub |>
    group_by(SNPS, MAPPED_TRAIT) |>
    slice_max(PVALUE_MLOG, n = 1, with_ties = FALSE) |>
    ungroup()

snp_trait <- sub[, c("SNPS", "MAPPED_TRAIT")]
snp_n <- sort(table(snp_trait$SNPS), decreasing = TRUE)
pleio_snps <- names(snp_n)[snp_n >= 2]

cat("Rows after dedup:", nrow(sub), "\n")
cat("Pleiotropic SNPs:", length(pleio_snps), "\n")

# Include all associations for the strongest pleiotropic SNPs, then fill the
# rest with a trait-balanced random sample. The documented flagship SNPs are
# always kept so the README examples remain valid.
top_pleio <- names(snp_n)[snp_n >= 2][seq_len(min(120, length(pleio_snps)))]
must_include <- c("rs814573", "rs429358")
top_pleio <- unique(c(top_pleio, must_include))

keep <- sub[sub$SNPS %in% top_pleio, ]

remaining <- sub[!(sub$SNPS %in% top_pleio), ]

target_size <- 2500
n_fill <- max(0, target_size - nrow(keep))

if (n_fill > 0 && nrow(remaining) > 0) {
    fill <- remaining |>
        group_by(MAPPED_TRAIT) |>
        sample_frac(min(1, n_fill / nrow(remaining))) |>
        ungroup()
    keep <- bind_rows(keep, fill)
}

gwas_subset <- keep |>
    arrange(CHR_ID, CHR_POS, SNPS, MAPPED_TRAIT) |>
    distinct(SNPS, MAPPED_TRAIT, .keep_all = TRUE)

gwas_subset <- as.data.table(gwas_subset)
rownames(gwas_subset) <- NULL

cat("\nFinal curated dataset\n")
cat("  rows:", nrow(gwas_subset), "\n")
cat("  SNPs:", length(unique(gwas_subset$SNPS)), "\n")
cat("  traits:", length(unique(gwas_subset$MAPPED_TRAIT)), "\n")

st <- table(unique(gwas_subset[, c("SNPS", "MAPPED_TRAIT")])$SNPS)
cat("  pleiotropic SNPs:", sum(st >= 2), "\n")
cat("  max traits per SNP:", max(st), "\n")
cat("  rs814573 present:", "rs814573" %in% gwas_subset$SNPS, "\n")
cat("  rs429358 present:", "rs429358" %in% gwas_subset$SNPS, "\n")

save(gwas_subset, file = "data/gwas_subset.rda", compress = "xz")
saveRDS(gwas_subset, file = "inst/extdata/gwas_subset.rds")

example <- head(as.data.frame(gwas_subset), 200)
write.table(example, file = "inst/extdata/example_gwas.tsv",
    sep = "\t", row.names = FALSE, quote = FALSE
)

cat("Wrote data/gwas_subset.rda, inst/extdata/gwas_subset.rds,",
    "inst/extdata/example_gwas.tsv\n"
)
