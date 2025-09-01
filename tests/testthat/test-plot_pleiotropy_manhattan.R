test_that("plot_pleiotropy_manhattan creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_pleiotropy_manhattan(pleio_data)
    expect_s3_class(p, "ggplot")
    expect_true("PLOT_POS" %in% names(p$data))
})

test_that("plot_pleiotropy_manhattan handles highlighting", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    if ("rs814573" %in% pleio_data$SNPS) {
        p <- plot_pleiotropy_manhattan(pleio_data, highlight_snp = "rs814573")
        expect_s3_class(p, "ggplot")
    }
})

test_that("plot_pleiotropy_manhattan validates input", {
    expect_error(
        plot_pleiotropy_manhattan("not_a_dataframe"),
        "Input must be a data.frame"
    )

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(
        plot_pleiotropy_manhattan(invalid_data),
        "Required columns missing"
    )

    empty_data <- data.frame(
        SNPS = character(), CHR_ID = character(),
        CHR_POS = character(), PVALUE_MLOG = numeric()
    )
    expect_error(plot_pleiotropy_manhattan(empty_data), "Input data is empty")
})
