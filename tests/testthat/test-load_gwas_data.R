test_that("load_gwas_data works correctly", {
    temp_file <- tempfile(fileext = ".tsv")
    test_data <- data.frame(
        SNPS = c("rs1", "rs2"),
        MAPPED_TRAIT = c("trait1", "trait2"),
        PVALUE_MLOG = c(8, 10)
    )
    write.table(test_data, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

    result <- load_gwas_data(temp_file)
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 2)
    expect_true(all(c("SNPS", "MAPPED_TRAIT", "PVALUE_MLOG") %in% names(result)))

    unlink(temp_file)
})

test_that("load_gwas_data handles missing files", {
    expect_error(load_gwas_data("nonexistent_file.txt"), "File not found")
    expect_error(
        load_gwas_data("../"),
        "Error reading file: File '../' is a directory. Not yet implemented."
    )
})

test_that("load_gwas_data handles different separators", {
    temp_file <- tempfile(fileext = ".csv")
    test_data <- data.frame(
        SNPS = c("rs1", "rs2"),
        MAPPED_TRAIT = c("trait1", "trait2")
    )
    write.table(test_data, temp_file, sep = ",", row.names = FALSE, quote = FALSE)

    result <- load_gwas_data(temp_file, sep = ",")
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 2)

    unlink(temp_file)
})
