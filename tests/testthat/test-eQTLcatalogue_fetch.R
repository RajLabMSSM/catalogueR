test_that("eQTLcatalogue_fetch works", {

    testthat::skip_on_cran()
    testthat::skip_if_offline()
    testthat::skip_if_not_installed("echodata")
    testthat::skip_if_not_installed("echotabix")

    ## Skip if eQTL Catalogue API is not reachable
    api_ok <- tryCatch({
        con <- url("https://www.ebi.ac.uk/eqtl/api/v3/", open = "r")
        close(con)
        TRUE
    }, error = function(e) FALSE)
    testthat::skip_if_not(api_ok,
                          message = "eQTL Catalogue API is not reachable")

    data("meta")
    dataset_id <- meta$unique_label[1]
    paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
    query_dat <- data.table::fread(paths$BST1)
    ## Detect position column name
    pos_col <- intersect(c("POS","BP"), colnames(query_dat))[1]
    testthat::skip_if(is.na(pos_col),
                      message = "No position column (POS/BP) found")

    #### Using data ####
    gwas.qtl <- catalogueR::eQTLcatalogue_fetch(
        dataset_id = dataset_id,
        query_granges = query_dat,
        method = "REST"
    )
    testthat::expect_true(methods::is(gwas.qtl, "data.table"))
    testthat::expect_gte(nrow(gwas.qtl), 1)

    #### using explicit range ####
    query_granges <- echotabix::construct_query(
      query_chrom = query_dat$CHR[[1]],
      query_start_pos = min(query_dat[[pos_col]]),
      query_end_pos = max(query_dat[[pos_col]]))
    GWAS.QTL_manual <- catalogueR::eQTLcatalogue_fetch(
      dataset_id = dataset_id,
      query_granges = query_granges)
    testthat::expect_true(methods::is(GWAS.QTL_manual, "data.table"))
    testthat::expect_gte(nrow(GWAS.QTL_manual), 1)
})
