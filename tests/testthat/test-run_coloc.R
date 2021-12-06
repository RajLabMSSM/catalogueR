test_that("run_coloc works", {
    
    gwas.qtl_paths <- example_eQTL_Catalogue_query_paths(
        fnames = c("BST1__Alasoo_2018.macrophage_IFNg+Salmonella.tsv.gz",
                   "BST1__Alasoo_2018.macrophage_IFNg.tsv.gz")
    )
    for(x in names(gwas.qtl_paths)){
        message(x)
        dat <- data.table::fread(gwas.qtl_paths[x], nThread = 1)
        testthat::expect_true(methods::is(dat,"data.table"))
        testthat::expect_gte(nrow(dat),30000)
    }
    
    coloc_QTLs <- catalogueR::run_coloc(gwas.qtl_paths = gwas.qtl_paths)

    testthat::expect_true(methods::is(coloc_QTLs, "data.table"))
    testthat::expect_gte(nrow(coloc_QTLs), 70000)
    testthat::expect_gte(nrow(subset(coloc_QTLs, PP.H4 > .1)), 9300)
})
