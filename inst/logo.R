library(hexSticker)
imgurl <-
    system.file("man/figures/pleior_logo02.png", package = "pleior")
sticker(
    imgurl,
    package = "",
    p_size = 12,
    s_x = 1,
    s_y = 1,
    s_width = .6,
    s_height = .8,
    p_color = "black",
    h_fill = "white",
    h_color = "white",
    filename = "man/figures/logo.png"
)
