test_that("eQTLcatalogue_query works", {
    
    library(dplyr)
    sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5)[seq_len(2)]
    qtl_search <- c("Alasoo_2018.macrophage_naive")

    #### Method 1: Merged results ####
    GWAS.QTL <- catalogueR::eQTLcatalogue_query(
        sumstats_paths = sumstats_paths,
        qtl_search = qtl_search,
        merge_with_gwas = TRUE,
        split_files = FALSE
    )
    testthat::expect_true(methods::is(GWAS.QTL, "data.table")) 
    testthat::expect_equal(dim(GWAS.QTL),c(80,56))
    testthat::expect_true(all(qtl_search %in% unique(GWAS.QTL$qtl_id)))

    #### Method 2: Split results ####
    gwas.qtl_paths <- catalogueR::eQTLcatalogue_query(
        sumstats_paths = sumstats_paths,
        qtl_search = qtl_search,
        merge_with_gwas = TRUE,
        split_files = TRUE
    )
    for (f in names(gwas.qtl_paths)) {
        print(f)
        dat <- data.table::fread(gwas.qtl_paths[[f]])
        testthat::expect_equal(dim(dat),c(55,53))
    }
    #### Gather files ####
    GWAS.QTL2 <- catalogueR::merge_files(file_paths = gwas.qtl_paths)
    testthat::expect_true(methods::is(GWAS.QTL2, "data.table"))
    testthat::expect_equal(nrow(GWAS.QTL), nrow(GWAS.QTL2))
})
