test_that("multiplication works", {
    expect_equal(2 * 2, 4)
})
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
