test_that("eQTLcatalogue_example_queries_local returns named paths", {
    result <- catalogueR:::eQTLcatalogue_example_queries_local()
    expect_true(is.character(result))
    # All results should be named
    expect_true(!is.null(names(result)))
})

test_that("eQTLcatalogue_example_queries_local with locus_names=TRUE", {
    result <- catalogueR:::eQTLcatalogue_example_queries_local(
        locus_names = TRUE
    )
    expect_true(is.character(result))
    # Names should be locus names (first part before underscore)
    if (length(result) > 0) {
        expect_true(!is.null(names(result)))
    }
})

test_that("eQTLcatalogue_example_queries_local with bad Rlib_path", {
    result <- catalogueR:::eQTLcatalogue_example_queries_local(
        Rlib_path = "/nonexistent/path"
    )
    # Should return empty character since the path doesn't exist
    expect_true(is.character(result))
    expect_length(result, 0)
})
