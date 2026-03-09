test_that("COLOC_heatmap works", {

    testthat::skip_on_cran()
    testthat::skip_if_offline()

    coloc_QTLs <- catalogueR::COLOC_get_example_res()
    gg_coloc <- catalogueR:: COLOC_heatmap(
        coloc_QTLs = coloc_QTLs,
        coloc_thresh = .5
    )
    testthat::expect_true(methods::is(gg_coloc$plot, "gg"))
})
