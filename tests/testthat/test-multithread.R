test_that("multithread_handler detects conflicts", {
    expect_warning(
        result <- catalogueR:::multithread_handler(
            multithread_qtl = TRUE,
            multithread_loci = TRUE,
            multithread_tabix = FALSE
        ),
        "Only one multithreading option"
    )
    expect_true(result$qtl)
    expect_false(result$loci)
    expect_false(result$tabix)
})

test_that("multithread_handler allows single option", {
    result <- catalogueR:::multithread_handler(
        multithread_qtl = TRUE,
        multithread_loci = FALSE,
        multithread_tabix = FALSE
    )
    expect_true(result$qtl)
    expect_false(result$loci)
    expect_false(result$tabix)
})

test_that("multithread_handler allows all FALSE", {
    result <- catalogueR:::multithread_handler(
        multithread_qtl = FALSE,
        multithread_loci = FALSE,
        multithread_tabix = FALSE
    )
    expect_false(result$qtl)
    expect_false(result$loci)
    expect_false(result$tabix)
})

test_that("multithread_optimizer chooses QTL when more datasets than loci", {
    result <- suppressMessages(
        catalogueR:::multithread_optimizer(
            qtl_datasets = c("ds1", "ds2", "ds3"),
            sumstats_paths = c("path1"),
            verbose = FALSE
        )
    )
    expect_true(result$qtl)
    expect_false(result$loci)
    expect_false(result$tabix)
})

test_that("multithread_optimizer chooses loci when more loci than datasets", {
    result <- suppressMessages(
        catalogueR:::multithread_optimizer(
            qtl_datasets = c("ds1"),
            sumstats_paths = c("path1", "path2", "path3"),
            verbose = FALSE
        )
    )
    expect_false(result$qtl)
    expect_true(result$loci)
    expect_false(result$tabix)
})
