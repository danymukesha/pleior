test_that("plot_pleiotropy_heatmap creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    p <- plot_pleiotropy_heatmap(pleio_data, top_n_snps = 3, top_n_traits = 3)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pleiotropy_heatmap validates input", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")

    expect_error(
        plot_pleiotropy_heatmap("not_a_dataframe"),
        "Input must be a data.frame"
    )

    invalid_data <- data.frame(wrong_col = 1:5)
    expect_error(
        plot_pleiotropy_heatmap(invalid_data),
        "Column .* not found in data"
    )
})

test_that("plot_regional_association creates valid plot", {
    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    skip_if(nrow(pleio_data) == 0, "No pleiotropic SNPs in test data")
    skip_if(!"rs814573" %in% pleio_data$SNPS, "Target SNP not in test data")

    p <- plot_regional_association(pleio_data, target_snp = "rs814573")
    expect_s3_class(p, "ggplot")
})

test_that("plot_regional_association validates input", {
    expect_error(
        plot_regional_association("not_a_dataframe", "rs814573"),
        "Input must be a data.frame"
    )

    expect_error(
        plot_regional_association(data.frame(SNPS = "rs1"), "rs814573"),
        "Required columns missing"
    )

    data(gwas_subset, envir = environment())
    pleio_data <- detect_pleiotropy(gwas_subset)

    expect_error(
        plot_regional_association(pleio_data),
        "target_snp is required"
    )
})
