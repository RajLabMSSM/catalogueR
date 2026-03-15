test_that("meta dataset is loadable and has expected structure", {
    data("meta", package = "catalogueR", envir = environment())
    expect_s3_class(meta, "data.frame")
    expect_true(nrow(meta) > 100)
    expected_cols <- c(
        "study_id", "dataset_id", "study_label",
        "sample_group", "tissue_label", "quant_method"
    )
    for (col in expected_cols) {
        expect_true(col %in% colnames(meta),
                    info = paste("Missing column:", col))
    }
})

test_that("eQTLcatalogue_header dataset is loadable", {
    data("eQTLcatalogue_header", package = "catalogueR", envir = environment())
    expect_true(is.character(eQTLcatalogue_header))
    expect_true(length(eQTLcatalogue_header) > 5)
})

test_that("meta contains Alasoo_2018 study", {
    data("meta", package = "catalogueR", envir = environment())
    expect_true("Alasoo_2018" %in% meta$study_label)
})

test_that("meta contains multiple quant methods", {
    data("meta", package = "catalogueR", envir = environment())
    expect_true(length(unique(meta$quant_method)) > 1)
})
