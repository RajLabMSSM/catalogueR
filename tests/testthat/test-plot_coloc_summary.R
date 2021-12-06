test_that("plot_coloc_summary works", {
    
    coloc_QTLs <- catalogueR::get_coloc_QTLs()
    gg_coloc <- catalogueR:: plot_coloc_summary(
        coloc_QTLs = coloc_QTLs,
        coloc_thresh = .5
    )
    testthat::expect_true(methods::is(gg_coloc, "gg"))
})
