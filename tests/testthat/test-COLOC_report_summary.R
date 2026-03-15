test_that("COLOC_report_summary processes results correctly", {
    testthat::skip_if_not_installed("coloc")

    # Create mock coloc results matching coloc::coloc.abf output structure
    coloc_res <- list(
        summary = c(
            nsnps = 100,
            PP.H0.abf = 0.01,
            PP.H1.abf = 0.02,
            PP.H2.abf = 0.03,
            PP.H3.abf = 0.04,
            PP.H4.abf = 0.90
        ),
        results = data.table::data.table(
            snp = c("rs1", "rs2", "rs3"),
            SNP.PP.H4 = c(0.6, 0.25, 0.05),
            pvalues.df1 = c(1e-5, 1e-3, 0.5),
            pvalues.df2 = c(1e-6, 1e-4, 0.1)
        )
    )
    result <- suppressMessages(
        catalogueR:::COLOC_report_summary(
            coloc.res = coloc_res,
            coloc_thresh = 0.8
        )
    )
    expect_s3_class(result, "data.table")
    expect_true("Colocalized" %in% colnames(result))
    expect_equal(nrow(result), 3)
})

test_that("COLOC_report_summary reports colocalization status", {
    testthat::skip_if_not_installed("coloc")

    # Colocalized case: PP.H3 + PP.H4 >= thresh and PP.H4/PP.H3 >= 2
    coloc_res_yes <- list(
        summary = c(
            nsnps = 50,
            PP.H0.abf = 0.01,
            PP.H1.abf = 0.01,
            PP.H2.abf = 0.01,
            PP.H3.abf = 0.07,
            PP.H4.abf = 0.90
        ),
        results = data.table::data.table(
            snp = "rs1",
            SNP.PP.H4 = 0.9,
            pvalues.df1 = 1e-5,
            pvalues.df2 = 1e-6
        )
    )
    expect_message(
        catalogueR:::COLOC_report_summary(
            coloc.res = coloc_res_yes,
            coloc_thresh = 0.8
        ),
        "colocalized"
    )
})

test_that("COLOC_report_summary handles NOT colocalized case", {
    testthat::skip_if_not_installed("coloc")

    coloc_res_no <- list(
        summary = c(
            nsnps = 50,
            PP.H0.abf = 0.50,
            PP.H1.abf = 0.20,
            PP.H2.abf = 0.20,
            PP.H3.abf = 0.05,
            PP.H4.abf = 0.05
        ),
        results = data.table::data.table(
            snp = "rs1",
            SNP.PP.H4 = 0.05,
            pvalues.df1 = 0.5,
            pvalues.df2 = 0.1
        )
    )
    expect_message(
        catalogueR:::COLOC_report_summary(
            coloc.res = coloc_res_no,
            coloc_thresh = 0.8
        ),
        "NOT colocalized"
    )
})
