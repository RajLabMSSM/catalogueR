test_that("COLOC_get_res runs coloc.abf on mock data", {
    testthat::skip_if_not_installed("coloc")

    set.seed(42)
    n_snps <- 50
    snps <- paste0("rs", seq_len(n_snps))
    mafs <- runif(n_snps, 0.05, 0.45)

    qtl_egene <- data.table::data.table(
        SNP = snps,
        pvalue.QTL = runif(n_snps, 1e-6, 0.9),
        beta.QTL = rnorm(n_snps, 0, 0.5),
        se.QTL = abs(rnorm(n_snps, 0.1, 0.05)),
        maf.QTL = mafs,
        an.QTL = rep(200, n_snps),
        gene.QTL = rep("TESTGENE", n_snps)
    )

    gwas_region <- data.table::data.table(
        SNP = snps,
        P = runif(n_snps, 1e-8, 0.9),
        Effect = rnorm(n_snps, 0, 0.3),
        StdErr = abs(rnorm(n_snps, 0.1, 0.05)),
        MAF = mafs,
        N = rep(10000, n_snps),
        Locus = rep("TestLocus", n_snps),
        CHR = rep(1, n_snps),
        POS = seq(100000, by = 1000, length.out = n_snps)
    )

    result <- suppressMessages(
        suppressWarnings(
            catalogueR:::COLOC_get_res(
                qtl.egene = qtl_egene,
                gwas.region = gwas_region,
                merge_by_rsid = TRUE,
                coloc_thresh = 0.8,
                method = "abf",
                verbose = FALSE
            )
        )
    )

    expect_true(is.list(result))
    expect_true("summary" %in% names(result))
    expect_true("results" %in% names(result))
    expect_true("Locus" %in% names(result))
    # summary should have PP values
    expect_true("PP.H0.abf" %in% names(result$summary))
    expect_true("PP.H4.abf" %in% names(result$summary))
    # results should be a data.frame with SNP-level PP
    expect_true(nrow(result$results) > 0)
})

test_that("COLOC_get_res handles no shared SNPs", {
    testthat::skip_if_not_installed("coloc")

    qtl_egene <- data.table::data.table(
        SNP = c("rs100", "rs200"),
        pvalue.QTL = c(0.01, 0.05),
        beta.QTL = c(0.5, 0.3),
        maf.QTL = c(0.2, 0.3),
        an.QTL = c(200, 200)
    )

    gwas_region <- data.table::data.table(
        SNP = c("rs999", "rs888"),
        P = c(1e-5, 0.5),
        Effect = c(0.3, 0.1),
        MAF = c(0.2, 0.3),
        N = c(10000, 10000),
        Locus = c("TestLocus", "TestLocus")
    )

    result <- suppressMessages(
        catalogueR:::COLOC_get_res(
            qtl.egene = qtl_egene,
            gwas.region = gwas_region,
            merge_by_rsid = TRUE,
            method = "abf",
            verbose = FALSE
        )
    )

    expect_true(is.list(result))
    # Should indicate no shared SNPs
    expect_true("summary" %in% names(result))
})

test_that("COLOC_get_res errors without N column", {
    testthat::skip_if_not_installed("coloc")

    qtl_egene <- data.table::data.table(
        SNP = "rs1",
        pvalue.QTL = 0.01,
        beta.QTL = 0.5,
        maf.QTL = 0.2,
        an.QTL = 200
    )
    # No N column in GWAS region
    gwas_region <- data.table::data.table(
        SNP = "rs1",
        P = 1e-5,
        Effect = 0.3,
        MAF = 0.2,
        Locus = "TestLocus"
    )

    expect_error(
        suppressMessages(
            catalogueR:::COLOC_get_res(
                qtl.egene = qtl_egene,
                gwas.region = gwas_region,
                verbose = FALSE
            )
        ),
        "N.*column.*not detected"
    )
})

test_that("COLOC_get_res removes duplicated SNPs", {
    testthat::skip_if_not_installed("coloc")

    set.seed(123)
    # Input with duplicates
    qtl_egene <- data.table::data.table(
        SNP = c("rs1", "rs1", "rs2"),
        pvalue.QTL = c(0.01, 0.02, 0.05),
        beta.QTL = c(0.5, 0.4, 0.3),
        se.QTL = c(0.1, 0.1, 0.1),
        maf.QTL = c(0.2, 0.2, 0.3),
        an.QTL = c(200, 200, 200)
    )

    gwas_region <- data.table::data.table(
        SNP = c("rs1", "rs2"),
        P = c(1e-5, 0.5),
        Effect = c(0.3, 0.1),
        StdErr = c(0.1, 0.1),
        MAF = c(0.2, 0.3),
        N = c(10000, 10000),
        Locus = c("TestLocus", "TestLocus"),
        CHR = c(1, 1),
        POS = c(100, 200)
    )

    result <- suppressMessages(
        suppressWarnings(
            catalogueR:::COLOC_get_res(
                qtl.egene = qtl_egene,
                gwas.region = gwas_region,
                method = "abf",
                verbose = FALSE
            )
        )
    )
    expect_true(is.list(result))
    # Should successfully run despite duplicate input
    expect_true("summary" %in% names(result))
})
