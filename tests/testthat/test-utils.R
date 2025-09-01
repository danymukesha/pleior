library(testthat)

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
})
