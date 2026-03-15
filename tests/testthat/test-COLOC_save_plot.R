test_that("COLOC_save_plot saves to correct path", {
    testthat::skip_if_not_installed("ggplot2")

    gg <- ggplot2::ggplot(
        data.frame(x = 1:5, y = 1:5),
        ggplot2::aes(x = x, y = y)
    ) + ggplot2::geom_point()

    coloc_dat <- data.table::data.table(
        dataset_id = c("ds1", "ds2"),
        Locus.eGene = c("BST1 (GENE1)", "LRRK2 (GENE2)")
    )

    tmp_dir <- tempfile("coloc_test_")
    dir.create(tmp_dir)

    result <- suppressMessages(
        suppressWarnings(
            catalogueR:::COLOC_save_plot(
                gg_coloc = gg,
                save_dir = tmp_dir,
                coloc_dat = coloc_dat,
                coloc_thresh = 0.8
            )
        )
    )
    expect_true(file.exists(result))
    expect_true(grepl("coloc_PP0.8.png", result))

    unlink(tmp_dir, recursive = TRUE)
})

test_that("COLOC_save_plot returns NULL when save_dir is NULL", {
    testthat::skip_if_not_installed("ggplot2")

    gg <- ggplot2::ggplot(
        data.frame(x = 1:5, y = 1:5),
        ggplot2::aes(x = x, y = y)
    ) + ggplot2::geom_point()

    coloc_dat <- data.table::data.table(
        dataset_id = "ds1",
        Locus.eGene = "BST1 (GENE1)"
    )

    result <- catalogueR:::COLOC_save_plot(
        gg_coloc = gg,
        save_dir = NULL,
        coloc_dat = coloc_dat,
        coloc_thresh = 0.8
    )
    expect_null(result)
})

test_that("COLOC_save_plot uses fixed dimensions at thresh=0.99", {
    testthat::skip_if_not_installed("ggplot2")

    gg <- ggplot2::ggplot(
        data.frame(x = 1:5, y = 1:5),
        ggplot2::aes(x = x, y = y)
    ) + ggplot2::geom_point()

    coloc_dat <- data.table::data.table(
        dataset_id = "ds1",
        Locus.eGene = "BST1 (GENE1)"
    )

    tmp_dir <- tempfile("coloc_test_99_")
    dir.create(tmp_dir)

    result <- suppressMessages(
        suppressWarnings(
            catalogueR:::COLOC_save_plot(
                gg_coloc = gg,
                save_dir = tmp_dir,
                coloc_dat = coloc_dat,
                coloc_thresh = 0.99
            )
        )
    )
    expect_true(file.exists(result))
    expect_true(grepl("coloc_PP0.99.png", result))

    unlink(tmp_dir, recursive = TRUE)
})
