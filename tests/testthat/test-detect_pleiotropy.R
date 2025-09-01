test_that("detect_pleiotropy identifies pleiotropic SNPs", {
    data(gwas_subset, envir = environment())

    result <- detect_pleiotropy(gwas_subset)
    expect_s3_class(result, "data.table")

    if (nrow(result) > 0) {
        expect_true(all(result$N_TRAITS > 1))
        expect_true("TRAITS" %in% names(result))
    }
})

test_that("detect_pleiotropy handles trait filtering", {
    data(gwas_subset, envir = environment())

    specific_traits <- c("Alzheimer disease", "myocardial infarction")
    result <- detect_pleiotropy(gwas_subset, traits = specific_traits)

    expect_s3_class(result, "data.table")
    if (nrow(result) > 0) {
        expect_true(all(grepl(
            paste(specific_traits, collapse = "|"),
            result$MAPPED_TRAIT
        )))
    }
})

test_that("detect_pleiotropy validates input", {
    expect_error(detect_pleiotropy("not_a_dataframe"), "Input must be a data.frame")

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(detect_pleiotropy(invalid_data), "Required columns missing")

    empty_data <- data.frame(
        SNPS = character(), MAPPED_TRAIT = character(),
        PVALUE_MLOG = numeric()
    )
    expect_error(detect_pleiotropy(empty_data), "Input data is empty")
})

test_that("detect_pleiotropy handles wrong trait", {
    test_data <- data.frame(
        SNPS = c("rs1", "rs2", "rs3"),
        MAPPED_TRAIT = c("trait1", "trait2", "trait3"),
        PVALUE_MLOG = c(8, NA, 10)
    )

    expect_error(
        detect_pleiotropy(test_data, traits = "wrong_trait"),
        "No data found for specified traits"
    )
    expect_warning(
        detect_pleiotropy(test_data),
        "No pleiotropic SNPs found with current parameters"
    )
})
