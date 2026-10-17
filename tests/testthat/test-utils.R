test_that("validate_gwas_data correctly validates GWAS data", {
    required_cols <- c("SNP", "P", "CHR", "BP")
    gwas_data_valid <- data.frame(
        SNP = "rs123", P = 0.05, CHR = 1, BP = 123456
    )
    expect_true(validate_gwas_data(gwas_data_valid, required_cols))
    expect_false(validate_gwas_data("not_a_df", required_cols))
    expect_false(validate_gwas_data(data.frame(), required_cols))
    gwas_data_missing_col <- data.frame(
        SNP = "rs123", P = 0.05, CHR = 1
    )
    expect_false(validate_gwas_data(gwas_data_missing_col, required_cols))
})

test_that("format_pvalues correctly formats p-values", {
    input <- c(1, 2, 3) # -log10(p) = 1 → p = 0.1, etc.
    expected <- formatC(10^(-input), format = "e", digits = 2)
    expect_equal(format_pvalues(input), expected)

    input <- 301 # -log10(p) = 301 → p < 1e-300
    expect_equal(format_pvalues(input), "< 1e-300")

    input <- c(1, 301)
    expected <- c(formatC(10^(-1), format = "e", digits = 2), "< 1e-300")
    expect_equal(format_pvalues(input), expected)

    expect_true(is.na(format_pvalues(NA)))
})

test_that("split_mapped_traits explodes multi-trait fields", {
    df <- data.frame(
        SNPS = c("rs1", "rs2", "rs3"),
        MAPPED_TRAIT = c("body mass index, obesity", "asthma", "a;b;c"),
        PVALUE_MLOG = c(8, 9, 10)
    )

    out <- split_mapped_traits(df)
    expect_equal(nrow(out), 6)
    expect_setequal(
        out$MAPPED_TRAIT,
        c("body mass index", "obesity", "asthma", "a", "b", "c")
    )
    expect_equal(out$PVALUE_MLOG[out$MAPPED_TRAIT == "obesity"], 8)
})

test_that("adjust_gwas_pvalues applies multiple-testing correction", {
    df <- data.frame(
        SNPS = paste0("rs", 1:5),
        MAPPED_TRAIT = "trait",
        PVALUE_MLOG = c(5, 6, 7, 8, 9)
    )

    expect_identical(adjust_gwas_pvalues(df, "none"), df)

    bonf <- adjust_gwas_pvalues(df, "bonferroni")
    expect_true(all(bonf$PVALUE_MLOG < df$PVALUE_MLOG))

    fdr <- adjust_gwas_pvalues(df, "fdr")
    expect_true(all(fdr$PVALUE_MLOG <= df$PVALUE_MLOG))
})

test_that("clump_by_position collapses co-located SNPs", {
    df <- data.frame(
        SNPS = c("rs1", "rs2", "rs3", "rs4"),
        MAPPED_TRAIT = "trait",
        PVALUE_MLOG = c(5, 9, 8, 6),
        CHR_ID = c("1", "1", "1", "2"),
        CHR_POS = c(100000, 100001, 5000000, 100000)
    )

    out <- clump_by_position(df, clump_kb = 10)
    # rs2 is the lead SNP for chr1 window near 100kb; rs1 is dropped
    expect_false("rs1" %in% out$SNPS)
    expect_true("rs2" %in% out$SNPS)
    # chromosome 2 SNP is retained independently
    expect_true("rs4" %in% out$SNPS)
})
