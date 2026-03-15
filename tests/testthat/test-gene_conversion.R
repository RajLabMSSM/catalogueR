test_that("ensembl_to_hgnc converts known IDs", {
    testthat::skip_if_not_installed("EnsDb.Hsapiens.v75")
    testthat::skip_if_not_installed("AnnotationDbi")

    result <- suppressMessages(
        catalogueR::ensembl_to_hgnc(
            ensembl_ids = c("ENSG00000176697", "ENSG00000128573"),
            verbose = FALSE
        )
    )
    expect_true(is.character(result))
    expect_length(result, 2)
    # ENSG00000176697 = BDNF, ENSG00000128573 = FOXP2
    expect_true(!all(is.na(result)))
})

test_that("ensembl_to_hgnc handles NA in input", {
    testthat::skip_if_not_installed("EnsDb.Hsapiens.v75")
    testthat::skip_if_not_installed("AnnotationDbi")

    result <- suppressMessages(
        catalogueR::ensembl_to_hgnc(
            ensembl_ids = c("ENSG00000176697", NA),
            verbose = FALSE
        )
    )
    expect_true(is.character(result))
})

test_that("ensembl_to_hgnc unique_only deduplicates", {
    testthat::skip_if_not_installed("EnsDb.Hsapiens.v75")
    testthat::skip_if_not_installed("AnnotationDbi")

    result <- suppressMessages(
        catalogueR::ensembl_to_hgnc(
            ensembl_ids = c("ENSG00000176697", "ENSG00000176697"),
            unique_only = TRUE,
            verbose = FALSE
        )
    )
    expect_length(result, 1)
})

test_that("hgnc_to_ensembl converts known symbols", {
    testthat::skip_if_not_installed("EnsDb.Hsapiens.v75")
    testthat::skip_if_not_installed("AnnotationDbi")

    result <- suppressMessages(
        catalogueR::hgnc_to_ensembl(
            gene_symbols = c("BDNF", "BST1"),
            verbose = FALSE
        )
    )
    expect_true(is.character(result))
    expect_length(result, 2)
    expect_true(all(grepl("^ENSG", result)))
})

test_that("hgnc_to_ensembl handles NA", {
    testthat::skip_if_not_installed("EnsDb.Hsapiens.v75")
    testthat::skip_if_not_installed("AnnotationDbi")

    result <- suppressMessages(
        catalogueR::hgnc_to_ensembl(
            gene_symbols = c("BDNF", NA),
            verbose = FALSE
        )
    )
    expect_true(is.character(result))
})
