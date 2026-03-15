test_that("eQTLcatalogue_get_header returns cached header", {
    result <- catalogueR:::eQTLcatalogue_get_header(force_new_header = FALSE)
    expect_true(is.character(result))
    expect_true(length(result) > 0)
    # Should contain known column names from eQTL Catalogue
    expect_true("pvalue" %in% result)
    expect_true("beta" %in% result)
    expect_true("rsid" %in% result)
    expect_true("gene_id" %in% result)
    expect_true("maf" %in% result)
})
