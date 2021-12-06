test_that("eQTL_Catalogue.fetch works", {
    
    data("meta")
    unique_id <- meta$unique_id[1]
    gwas_data <- echodata::BST1[seq_len(2),]
    
    #### Using data ####
    gwas.qtl <- catalogueR::eQTL_Catalogue.fetch(
        unique_id = unique_id,
        gwas_data = gwas_data, 
        method = "REST"
    )
    testthat::expect_true(methods::is(gwas.qtl, "data.table"))
    testthat::expect_gte(nrow(gwas.qtl), 1000)
    
    #### using explicit range ####
    GWAS.QTL_manual <- catalogueR::eQTL_Catalogue.fetch(
        unique_id = unique_id,
        chrom = gwas_data$CHR[1],
        bp_lower = min(gwas_data$POS),
        bp_upper = max(gwas_data$POS))
    testthat::expect_true(methods::is(GWAS.QTL_manual, "data.table"))
    testthat::expect_gte(nrow(GWAS.QTL_manual), 1000)
    
    #### Compare across methods ####
    testthat::expect_true(all.equal(gwas.qtl, 
                                    GWAS.QTL_manual))
})
