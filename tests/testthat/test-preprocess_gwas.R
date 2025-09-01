test_that("preprocess_gwas filters correctly", {
    data(gwas_subset, envir = environment())

    result <- preprocess_gwas(gwas_subset, pvalue_threshold = 1e-5)
    expect_s3_class(result, "data.table")
    expect_true(all(result$PVALUE_MLOG >= -log10(1e-5)))
})

test_that("preprocess_gwas handles missing data", {
    test_data <- data.frame(
        SNPS = c("rs1", "rs2", "rs3"),
        MAPPED_TRAIT = c("trait1", "trait2", "trait3"),
        PVALUE_MLOG = c(8, NA, 10)
    )

    result <- preprocess_gwas(test_data)
    expect_equal(nrow(result), 2)
    expect_false(any(is.na(result$PVALUE_MLOG)))
})

test_that("preprocess_gwas validates input", {
    expect_error(preprocess_gwas("not_a_dataframe"), "Input must be a data.frame")
    expect_error(preprocess_gwas(data.frame()), "Input data is empty")

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(preprocess_gwas(invalid_data), "Required columns missing")
})

test_that("preprocess_gwas handles wrong trait", {
    test_data <- data.frame(
        SNPS = c("rs1", "rs2", "rs3"),
        MAPPED_TRAIT = c("trait1", "trait2", "trait3"),
        PVALUE_MLOG = c(8, NA, 10)
    )

    expect_error(
        preprocess_gwas(test_data, columns = "wrong_column"),
        "None of the specified columns found in data"
    )
    expect_warning(
        preprocess_gwas(test_data, pvalue_threshold = 0)
    )
})
