test_that("construct_locus_name builds name from POS column", {
    dat <- data.table::data.table(
        CHR = c(4, 4, 4),
        POS = c(100, 200, 300),
        SNP = paste0("rs", 1:3)
    )
    result <- suppressMessages(
        catalogueR:::construct_locus_name(dat, verbose = FALSE)
    )
    expect_equal(result, "locus_chr4-100-300")
})

test_that("construct_locus_name works with BP column", {
    dat <- data.table::data.table(
        CHR = c(1, 1),
        BP = c(100, 200),
        SNP = paste0("rs", 1:2)
    )
    result <- suppressMessages(
        catalogueR:::construct_locus_name(dat, verbose = FALSE)
    )
    expect_equal(result, "locus_chr1-100-200")
})

test_that("construct_locus_name errors without position column", {
    dat <- data.table::data.table(
        CHR = c(1, 1),
        SNP = paste0("rs", 1:2)
    )
    expect_error(
        catalogueR:::construct_locus_name(dat, verbose = FALSE),
        "No position column"
    )
})

test_that("parse_gwas.qtl_path extracts dataset_id", {
    path <- "/some/dir/BST1__Alasoo_2018.macrophage_naive.tsv.gz"
    result <- catalogueR:::parse_gwas.qtl_path(
        gwas.qtl_path = path,
        get_locus = FALSE,
        get_dataset_id = TRUE
    )
    expect_equal(result, "Alasoo_2018.macrophage_naive")
})

test_that("parse_gwas.qtl_path extracts locus", {
    path <- "/some/dir/LRRK2__Alasoo_2018.macrophage_IFNg.tsv.gz"
    result <- catalogueR:::parse_gwas.qtl_path(
        gwas.qtl_path = path,
        get_locus = TRUE,
        get_dataset_id = FALSE
    )
    expect_equal(result, "LRRK2")
})

test_that("parse_gwas.qtl_path extracts both", {
    path <- "MEX3C__GTEx_V8.Brain_Cortex.tsv.gz"
    result <- catalogueR:::parse_gwas.qtl_path(
        gwas.qtl_path = path,
        get_locus = TRUE,
        get_dataset_id = TRUE
    )
    expect_true(is.list(result))
    expect_equal(result$locus, "MEX3C")
    expect_equal(result$dataset_id, "GTEx_V8.Brain_Cortex")
})

test_that("make_split_path creates correct path structure", {
    tmp <- tempdir()
    result <- catalogueR:::make_split_path(
        output_dir = tmp,
        dataset_id = "Alasoo_2018.macrophage_naive",
        loc = "BST1"
    )
    expect_true(endsWith(result, "BST1__Alasoo_2018.macrophage_naive.tsv.gz"))
    expect_true(grepl("Alasoo_2018.macrophage_naive", dirname(result)))
    expected_name <- "BST1__Alasoo_2018.macrophage_naive"
    expect_equal(names(result), expected_name)
})

test_that("reconstruct_file_names builds correct file names", {
    coloc_QTLs <- data.table::data.table(
        Locus.GWAS = c("BST1", "BST1", "LRRK2"),
        dataset_id = c("ds1", "ds1", "ds2"),
        PP.H4 = c(0.9, 0.8, 0.7)
    )
    result <- catalogueR:::reconstruct_file_names(coloc_QTLs)
    expect_length(result, 2)
    expect_true("BST1__ds1.tsv.gz" %in% result)
    expect_true("LRRK2__ds2.tsv.gz" %in% result)
})
