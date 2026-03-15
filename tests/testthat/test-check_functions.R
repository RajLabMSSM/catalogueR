test_that("check_coord_input raises error when no coords provided", {
    expect_error(
        catalogueR:::check_coord_input(
            query_dat = NULL,
            chrom = NULL,
            bp_lower = NULL,
            bp_upper = NULL
        ),
        "User must specify coordinates"
    )
})

test_that("check_coord_input does not error when query_dat provided", {
    dat <- data.table::data.table(CHR = 1, POS = 100)
    expect_no_error(
        catalogueR:::check_coord_input(
            query_dat = dat,
            chrom = NULL,
            bp_lower = NULL,
            bp_upper = NULL
        )
    )
})

test_that("check_coord_input does not error when chrom/bp provided", {
    expect_no_error(
        catalogueR:::check_coord_input(
            query_dat = NULL,
            chrom = 1,
            bp_lower = 100,
            bp_upper = 200
        )
    )
})

test_that("check_dim prints dimensions", {
    df <- data.table::data.table(a = 1:10, b = letters[1:10])
    expect_message(
        catalogueR:::check_dim(df),
        "Data dimensions:"
    )
})

test_that("check_dim handles NULL gracefully", {
    # Should not error - has tryCatch
    expect_no_error(
        catalogueR:::check_dim(NULL)
    )
})

test_that("check_maf filters invalid MAF values", {
    gwas <- data.table::data.table(
        SNP = paste0("rs", 1:5),
        MAF = c(0.1, 0.0, 0.5, NA, 1.0)
    )
    eqtl <- data.table::data.table(
        SNP = paste0("rs", 1:5),
        maf.QTL = c(0.2, 0.3, 0.0, 0.4, NA)
    )
    result <- catalogueR:::check_maf(gwas, eqtl, verbose = FALSE)

    # MAF=0, MAF=1, MAF=NA should be removed
    expect_true(all(result$gwas_shared$MAF > 0))
    expect_true(all(result$gwas_shared$MAF < 1))
    expect_true(!any(is.na(result$gwas_shared$MAF)))

    expect_true(all(result$eqtl_shared$maf.QTL > 0))
    expect_true(all(result$eqtl_shared$maf.QTL < 1))
    expect_true(!any(is.na(result$eqtl_shared$maf.QTL)))
})

test_that("check_maf borrows MAF from QTL when GWAS missing", {
    gwas <- data.table::data.table(
        SNP = paste0("rs", 1:3),
        P = c(0.01, 0.05, 0.001)
        # No MAF column
    )
    eqtl <- data.table::data.table(
        SNP = paste0("rs", 1:3),
        maf.QTL = c(0.1, 0.2, 0.3)
    )
    expect_message(
        result <- catalogueR:::check_maf(gwas, eqtl, verbose = TRUE),
        "Borrowing MAF from QTL"
    )
    expect_true("MAF" %in% colnames(result$gwas_shared))
})

test_that("check_maf converts MAF to numeric", {
    gwas <- data.table::data.table(
        SNP = paste0("rs", 1:2),
        MAF = c("0.1", "0.3")
    )
    eqtl <- data.table::data.table(
        SNP = paste0("rs", 1:2),
        maf.QTL = c("0.2", "0.4")
    )
    result <- catalogueR:::check_maf(gwas, eqtl, verbose = FALSE)
    expect_true(is.numeric(result$gwas_shared$MAF))
    expect_true(is.numeric(result$eqtl_shared$maf.QTL))
})

test_that("check_sample_size adds N column when missing", {
    dat <- data.table::data.table(SNP = "rs1", P = 0.05)
    result <- catalogueR:::check_sample_size(dat, sample_size = 1000)
    expect_true("N" %in% colnames(result))
    expect_equal(result$N, 1000)
})

test_that("check_sample_size keeps existing N column", {
    dat <- data.table::data.table(SNP = "rs1", P = 0.05, N = 500)
    result <- catalogueR:::check_sample_size(dat, sample_size = 1000)
    expect_equal(result$N, 500)
})

test_that("check_sample_size errors when N missing and no sample_size", {
    dat <- data.table::data.table(SNP = "rs1", P = 0.05)
    expect_error(
        catalogueR:::check_sample_size(dat, sample_size = NULL),
        "N.*column.*sample size.*not detected"
    )
})

test_that("check_missing_queries identifies missing queries", {
    GWAS.QTL <- data.table::data.table(
        Locus.GWAS = c("BST1", "BST1"),
        dataset_id = c("Alasoo_2018.macrophage_naive",
                        "Alasoo_2018.macrophage_naive")
    )
    missing <- suppressMessages(
        catalogueR:::check_missing_queries(
            qtl_datasets = c("Alasoo_2018.macrophage_naive",
                             "Alasoo_2018.macrophage_IFNg"),
            loci = c("BST1", "LRRK2"),
            GWAS.QTL = GWAS.QTL,
            verbose = FALSE
        )
    )
    # BST1__Alasoo_2018.macrophage_IFNg, LRRK2__* should be missing
    expect_true(length(missing) > 0)
    expect_true("BST1__Alasoo_2018.macrophage_IFNg" %in% missing)
    expect_true("LRRK2__Alasoo_2018.macrophage_naive" %in% missing)
    expect_true("LRRK2__Alasoo_2018.macrophage_IFNg" %in% missing)
})

test_that("check_missing_queries returns empty when all found", {
    GWAS.QTL <- data.table::data.table(
        Locus.GWAS = c("BST1", "BST1"),
        dataset_id = c("ds1", "ds1")
    )
    missing <- suppressMessages(
        catalogueR:::check_missing_queries(
            qtl_datasets = "ds1",
            loci = "BST1",
            GWAS.QTL = GWAS.QTL,
            verbose = FALSE
        )
    )
    expect_length(missing, 0)
})
