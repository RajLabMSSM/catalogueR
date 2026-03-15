test_that("merge_files merges CSV files into one data.table", {
    tmp1 <- tempfile(fileext = ".csv")
    tmp2 <- tempfile(fileext = ".csv")
    data.table::fwrite(data.table::data.table(x = 1:3, y = letters[1:3]), tmp1)
    data.table::fwrite(data.table::data.table(x = 4:6, y = letters[4:6]), tmp2)

    result <- suppressMessages(
        catalogueR::merge_files(
            file_paths = c(tmp1, tmp2),
            nThread = 1,
            verbose = FALSE
        )
    )
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 6)
    expect_true("file" %in% colnames(result))
    expect_true("x" %in% colnames(result))

    unlink(c(tmp1, tmp2))
})

test_that("merge_files handles RDS files", {
    tmp <- tempfile(fileext = ".rds")
    saveRDS(data.table::data.table(a = 1:2, b = 3:4), tmp)

    result <- suppressMessages(
        catalogueR::merge_files(
            file_paths = tmp,
            nThread = 1,
            verbose = FALSE
        )
    )
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 2)
    expect_true("file" %in% colnames(result))

    unlink(tmp)
})

test_that("merge_files with fill=TRUE handles mismatched columns", {
    tmp1 <- tempfile(fileext = ".csv")
    tmp2 <- tempfile(fileext = ".csv")
    data.table::fwrite(data.table::data.table(x = 1:2, y = "a"), tmp1)
    data.table::fwrite(data.table::data.table(x = 3:4, z = "b"), tmp2)

    result <- suppressMessages(
        catalogueR::merge_files(
            file_paths = c(tmp1, tmp2),
            nThread = 1,
            verbose = FALSE
        )
    )
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 4)
    # Both y and z should exist with NAs filled
    expect_true("y" %in% colnames(result))
    expect_true("z" %in% colnames(result))

    unlink(c(tmp1, tmp2))
})

test_that("merge_gwas_qtl merges on SNP column", {
    gwas <- data.table::data.table(
        SNP = c("rs1", "rs2", "rs3"),
        P = c(1e-8, 0.05, 0.01),
        CHR = c(1, 1, 1),
        POS = c(100, 200, 300)
    )
    qtl <- data.table::data.table(
        rsid.QTL = c("rs1", "rs2", "rs4"),
        pvalue.QTL = c(1e-5, 0.1, 0.02),
        beta.QTL = c(0.5, -0.2, 0.3)
    )
    result <- suppressMessages(
        catalogueR:::merge_gwas_qtl(
            query_dat = gwas,
            qtl.subset = qtl,
            verbose = FALSE
        )
    )
    expect_s3_class(result, "data.table")
    # Only rs1 and rs2 should match
    expect_equal(nrow(result), 2)
    expect_true(all(c("SNP", "P", "pvalue.QTL") %in% colnames(result)))
})

test_that("merge_gwas_qtl adds allele flip columns when A1/A2 present", {
    gwas <- data.table::data.table(
        SNP = c("rs1", "rs2"),
        P = c(1e-8, 0.05),
        A1 = c("A", "T"),
        A2 = c("G", "C")
    )
    qtl <- data.table::data.table(
        rsid.QTL = c("rs1", "rs2"),
        pvalue.QTL = c(1e-5, 0.1),
        ref.QTL = c("A", "G"),
        alt.QTL = c("G", "C")
    )
    result <- suppressMessages(
        catalogueR:::merge_gwas_qtl(
            query_dat = gwas,
            qtl.subset = qtl,
            verbose = FALSE
        )
    )
    expect_true("effect.is.ref" %in% colnames(result))
    expect_true("effect.is.alt" %in% colnames(result))
    # rs1: A1=A, ref=A => effect.is.ref=TRUE
    expect_true(result[SNP == "rs1"]$effect.is.ref)
})

test_that("merge_gwas_qtl returns qtl.subset on merge error", {
    # Bad data that would fail merge
    gwas <- "not a data.table"
    qtl <- data.table::data.table(
        rsid.QTL = "rs1",
        pvalue.QTL = 0.05
    )
    result <- suppressMessages(
        catalogueR:::merge_gwas_qtl(
            query_dat = gwas,
            qtl.subset = qtl,
            verbose = FALSE
        )
    )
    # Should fall back to returning qtl.subset
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 1)
})
