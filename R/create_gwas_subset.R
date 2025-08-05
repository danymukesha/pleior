library(data.table)

gwas_data <- data.table(
    SNPS = c(
        "rs814573", "rs814573", "rs7412", "rs7412", "rs2817462", "rs2817462",
        "rs10182181", "rs10182181", "rs3739081", "rs3739081"
    ),
    MAPPED_TRAIT = c(
        "Alzheimer disease", "myocardial infarction",
        "Alzheimer disease", "LDL cholesterol",
        "memory performance", "tyrosine measurement",
        "body mass index", "tyrosine measurement",
        "body mass index", "C-reactive protein"
    ),
    PVALUE_MLOG = c(
        672.699, 15.000, 122.3979, 9629, 6.30103, 9.69897,
        29.69897, 11.69897, 8.69897, 11.000
    ),
    CHR_ID = c("19", "19", "19", "19", "6", "6", "2", "2", "2", "2"),
    CHR_POS = c(
        "44908684", "44908684", "44919689", "44919689", "156588550",
        "156588550", "24927427", "24927427", "26732753", "26732753"
    )
)

gwas_subset <- gwas_data
usethis::use_data(gwas_subset, overwrite = TRUE)
