make_venn_data <- function() {
    data.frame(
        SNPS = c(
            "rs1", "rs1", "rs2", "rs2", "rs3", "rs3", "rs4", "rs4", "rs5", "rs5"
        ),
        MAPPED_TRAIT = c(
            "t1", "t2", "t1", "t3", "t2", "t4", "t3", "t4", "t1", "t4"
        ),
        PVALUE_MLOG = rep(10, 10),
        CHR_ID = rep("1", 10),
        CHR_POS = seq(1, 10)
    )
}

test_that("plot_venn_diagram supports 2, 3 and 4 sets", {
    df <- make_venn_data()

    p2 <- plot_venn_diagram(df, traits = c("t1", "t2"))
    expect_s3_class(p2, "ggplot")

    p3 <- plot_venn_diagram(df, traits = c("t1", "t2", "t3"))
    expect_s3_class(p3, "ggplot")

    p4 <- plot_venn_diagram(df, traits = c("t1", "t2", "t3", "t4"))
    expect_s3_class(p4, "ggplot")
    # the four sets are drawn as ellipse/circle layers
    expect_true(any(vapply(
        p4$layers, function(l) inherits(l$geom, "GeomCircle"), logical(1)
    )))
})

test_that("plot_venn_diagram validates input", {
    df <- make_venn_data()

    expect_error(
        plot_venn_diagram("not_a_dataframe", traits = c("t1", "t2")),
        "Input must be a data.frame"
    )
    expect_error(
        plot_venn_diagram(df),
        "Provide 2-4 traits"
    )
    expect_error(
        plot_venn_diagram(df, traits = "t1"),
        "Provide 2-4 traits"
    )
    expect_error(
        plot_venn_diagram(df, traits = c("t1", "t2", "t3", "t4", "t5")),
        "Provide 2-4 traits"
    )

    empty <- df[0, ]
    expect_error(
        plot_venn_diagram(empty, traits = c("t1", "t2")),
        "Input data is empty"
    )
})

test_that("plot_pleiotropy_landscape creates valid plot", {
    df <- make_venn_data()

    p <- plot_pleiotropy_landscape(df)
    expect_s3_class(p, "ggplot")

    expect_error(
        plot_pleiotropy_landscape("not_a_dataframe"),
        "Input must be a data.frame"
    )
    expect_error(
        plot_pleiotropy_landscape(data.frame(x = 1)),
        "Column 'PVALUE_MLOG' not found in data"
    )
})
