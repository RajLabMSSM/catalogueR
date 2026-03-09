test_that("eQTLcatalogue_query works", {

    testthat::skip_on_cran()
    testthat::skip_if_offline()
    testthat::skip_if_not_installed("echodata")
    testthat::skip_if_not_installed("dplyr")

    ## Skip if eQTL Catalogue API is not reachable
    api_ok <- tryCatch({
        con <- url("https://www.ebi.ac.uk/eqtl/api/v3/", open = "r")
        close(con)
        TRUE
    }, error = function(e) FALSE)
    testthat::skip_if_not(api_ok,
                          message = "eQTL Catalogue API is not reachable")

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
    testthat::expect_gte(nrow(GWAS.QTL), 0)
    if(nrow(GWAS.QTL) > 0) {
        testthat::expect_true(all(qtl_search %in% unique(GWAS.QTL$dataset_id)))
    }

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
        testthat::expect_true(methods::is(dat, "data.table"))
    }
    #### Gather files ####
    GWAS.QTL2 <- catalogueR::merge_files(file_paths = gwas.qtl_paths)
    testthat::expect_true(methods::is(GWAS.QTL2, "data.table"))
})
