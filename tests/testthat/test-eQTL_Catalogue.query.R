test_that("eQTL_Catalogue.query works", {
    
    library(dplyr)
    sumstats_paths <- catalogueR::example_sumstats_paths()[seq_len(2)]
    qtl_search <- c("Alasoo_2018.macrophage_naive")

    #### Method 1: Merged results ####
    GWAS.QTL <- catalogueR:: eQTL_Catalogue.query(
        sumstats_paths = sumstats_paths,
        qtl_search = qtl_search,
        merge_with_gwas = TRUE,
        split_files = FALSE
    )
    testthat::expect_true(methods::is(GWAS.QTL, "data.table"))
    counts <- GWAS.QTL %>%
        dplyr::group_by(qtl_id) %>%
        dplyr::count()
    testthat::expect_true(all(counts$n > 50000))
    testthat::expect_true(all(qtl_search
    %in% unique(GWAS.QTL$qtl_id)))

    #### Method 2: Split results ####
    gwas.qtl_paths <- catalogueR:: eQTL_Catalogue.query(
        sumstats_paths = sumstats_paths,
        qtl_search = qtl_search,
        merge_with_gwas = TRUE,
        split_files = TRUE
    )
    for (f in names(gwas.qtl_paths)) {
        print(f)
        dat <- data.table::fread(gwas.qtl_paths[[f]], nThread = 1)
        testthat::expect_gte(nrow(dat), 2000)
    }
    #### Gather files ####
    GWAS.QTL2 <- catalogueR::merge_files(file_paths = gwas.qtl_paths)
    testthat::expect_true(methods::is(GWAS.QTL2, "data.table"))
    testthat::expect_equal(nrow(GWAS.QTL), nrow(GWAS.QTL2))
})
