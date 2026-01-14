test_that("plot_pleiotropy_network creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_pleiotropy_network(pleio_data, top_n_snps = 3)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pleiotropy_network validates input", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    expect_error(
        plot_pleiotropy_network("not_a_dataframe"),
        "Input must be a data.frame"
    )

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(
        plot_pleiotropy_network(invalid_data),
        "Required columns missing"
    )
})

test_that("plot_trait_cooccurrence_network creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_trait_cooccurrence_network(pleio_data, top_n_traits = 3)
    expect_s3_class(p, "ggplot")
})

test_that("plot_trait_cooccurrence_network validates input", {
    expect_error(
        plot_trait_cooccurrence_network("not_a_dataframe"),
        "Input must be a data.frame"
    )

    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    if (nrow(pleio_data) > 0) {
        expect_error(
            plot_trait_cooccurrence_network(pleio_data, top_n_traits = 0),
            "Need at least 2 traits"
        )
    }
})
