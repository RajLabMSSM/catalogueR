test_that("COLOC_heatmap works", {
    
    coloc_QTLs <- catalogueR::COLOC_get_example_res()
    gg_coloc <- catalogueR:: COLOC_heatmap(
        coloc_QTLs = coloc_QTLs,
        coloc_thresh = .5
    )
    testthat::expect_true(methods::is(gg_coloc, "gg"))
})
