test_that("plot_multi_trait_manhattan creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_multi_trait_manhattan(pleio_data, max_traits = 2)
    expect_s3_class(p, "ggplot")
})

test_that("plot_multi_trait_manhattan validates input", {
    expect_error(
        plot_multi_trait_manhattan("not_a_dataframe"),
        "Input must be a data.frame"
    )

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(
        plot_multi_trait_manhattan(invalid_data),
        "Required columns missing"
    )
})

test_that("plot_effect_size_comparison creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_effect_size_comparison(pleio_data, top_n_snps = 3)
    expect_s3_class(p, "ggplot")
})

test_that("plot_effect_size_comparison validates input", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    expect_error(
        plot_effect_size_comparison("not_a_dataframe"),
        "Input must be a data.frame"
    )

    expect_error(
        plot_effect_size_comparison(pleio_data, effect_col = "nonexistent"),
        "Column .* not found in data"
    )
})

test_that("plot_effect_size_comparison honours color_by and error bars", {
    test_data <- data.frame(
        SNPS = rep(c("rs1", "rs2"), each = 2),
        MAPPED_TRAIT = rep(c("t1", "t2"), 2),
        PVALUE_MLOG = c(10, 11, 12, 13),
        CHR_ID = c("1", "1", "2", "2"),
        SE = c(0.1, 0.2, 0.15, 0.25)
    )

    expect_s3_class(
        plot_effect_size_comparison(test_data, color_by = "trait"),
        "ggplot"
    )
    expect_s3_class(
        plot_effect_size_comparison(test_data, color_by = "chromosome"),
        "ggplot"
    )
    expect_s3_class(
        plot_effect_size_comparison(test_data,
            show_error_bars = TRUE, error_col = "SE"
        ),
        "ggplot"
    )

    expect_error(
        plot_effect_size_comparison(test_data, show_error_bars = TRUE),
        "error_col must be specified"
    )
})
