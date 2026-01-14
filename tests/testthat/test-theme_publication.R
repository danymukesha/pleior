test_that("theme_pleiotropy_publication creates valid theme", {
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
        ggplot2::geom_point() +
        theme_pleiotropy_publication()

    expect_s3_class(p, "ggplot")
})

test_that("theme_pleiotropy_publication supports journal styles", {
    journal_styles <- c("default", "nature", "science", "pnas")

    for (style in journal_styles) {
        p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
            ggplot2::geom_point() +
            theme_pleiotropy_publication(journal_style = style)
        expect_s3_class(p, "ggplot")
    }
})

test_that("get_pleiotropy_colors returns valid colors", {
    colors_okabe_ito <- get_pleiotropy_colors("okabe_ito")
    expect_type(colors_okabe_ito, "character")
    expect_length(colors_okabe_ito, 10)

    colors_viridis <- get_pleiotropy_colors("viridis", n = 5)
    expect_type(colors_viridis, "character")
    expect_length(colors_viridis, 5)
})

test_that("get_pleiotropy_colors validates input", {
    expect_error(
        get_pleiotropy_colors("invalid_palette"),
        "Invalid palette_name"
    )
})
