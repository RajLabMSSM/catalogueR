test_that("eQTLcatalogue_fetch works", {
    
    data("meta")
    unique_id <- meta$unique_id[1]
    paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
    query_dat <- data.table::fread(paths$BST1)
    
    #### Using data ####
    gwas.qtl <- catalogueR::eQTLcatalogue_fetch(
        unique_id = unique_id,
        query_granges = query_dat, 
        method = "REST"
    )
    testthat::expect_true(methods::is(gwas.qtl, "data.table"))
    testthat::expect_gte(nrow(gwas.qtl), 1200)
    
    #### using explicit range ####
    query_granges <- echotabix::construct_query(
      query_chrom = query_dat$CHR[[1]],
      query_start_pos = min(query_dat$POS), 
      query_end_pos = max(query_dat$POS))
    GWAS.QTL_manual <- catalogueR::eQTLcatalogue_fetch(
      unique_id = unique_id, 
      query_granges = query_granges)
    testthat::expect_true(methods::is(GWAS.QTL_manual, "data.table"))
    testthat::expect_gte(nrow(GWAS.QTL_manual), 1200)
    
    #### Compare across methods ####
    testthat::expect_true(all.equal(gwas.qtl, 
                                    GWAS.QTL_manual))
    
    
    #### Using tabix #### 
    ## WARNING: if you test this iteratively, will eventually get blocked 
    ## by the EMBL-EBI server.
    gwas.qtl <- catalogueR::eQTLcatalogue_fetch(
      unique_id = unique_id,
      query_granges = query_dat, 
      method = "tabix"
    )
    testthat::expect_true(methods::is(gwas.qtl, "data.table"))
    testthat::expect_gte(nrow(gwas.qtl), 1200)
})
