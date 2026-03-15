test_that("eQTLcatalogue_list_datasets returns meta data", {
    meta <- suppressMessages(
        catalogueR::eQTLcatalogue_list_datasets(
            force_new = FALSE,
            verbose = FALSE
        )
    )
    expect_s3_class(meta, "data.frame")
    expect_true(nrow(meta) > 0)
    expect_true("dataset_id" %in% colnames(meta))
    expect_true("study_label" %in% colnames(meta))
    expect_true("tissue_label" %in% colnames(meta))
    expect_true("quant_method" %in% colnames(meta))
})

test_that("eQTLcatalogue_list_datasets adds System and Tissue_group", {
    meta <- suppressMessages(
        catalogueR::eQTLcatalogue_list_datasets(
            force_new = FALSE,
            verbose = FALSE
        )
    )
    expect_true("System" %in% colnames(meta))
    expect_true("Tissue_group" %in% colnames(meta))
    # Blood tissues should be categorized as "Blood"
    blood_rows <- meta[meta$tissue_label == "macrophage", ]
    if (nrow(blood_rows) > 0) {
        expect_true(all(blood_rows$System == "Blood"))
    }
})

test_that("eQTLcatalogue_list_datasets creates unique_label", {
    meta <- suppressMessages(
        catalogueR::eQTLcatalogue_list_datasets(
            force_new = FALSE,
            verbose = FALSE
        )
    )
    expect_true("unique_label" %in% colnames(meta))
    # unique_label should be study_label.sample_group
    first <- meta[1, ]
    expect_equal(
        first$unique_label,
        paste0(first$study_label, ".", first$sample_group)
    )
})

test_that("eQTLcatalogue_search_metadata returns all with NULL search", {
    result <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = NULL,
            verbose = FALSE
        )
    )
    expect_true(length(result) > 10)
})

test_that("eQTLcatalogue_search_metadata filters by search term", {
    result <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "Alasoo_2018",
            verbose = FALSE
        )
    )
    expect_true(length(result) > 0)
    expect_true(length(result) < 50)  # Should be a subset
})

test_that("eQTLcatalogue_search_metadata handles multiple search terms", {
    result_single <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "Alasoo_2018",
            verbose = FALSE
        )
    )
    result_multi <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = c("Alasoo_2018", "BLUEPRINT"),
            verbose = FALSE
        )
    )
    expect_true(length(result_multi) >= length(result_single))
})

test_that("eQTLcatalogue_search_metadata is case-insensitive", {
    result_lower <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "alasoo_2018",
            return_dataset_ids = FALSE,
            verbose = FALSE
        )
    )
    result_mixed <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "Alasoo_2018",
            return_dataset_ids = FALSE,
            verbose = FALSE
        )
    )
    expect_equal(length(result_lower), length(result_mixed))
})

test_that("eQTLcatalogue_search_metadata return_dataset_ids works", {
    result_ids <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "Alasoo_2018",
            return_dataset_ids = TRUE,
            verbose = FALSE
        )
    )
    result_labels <- suppressMessages(
        catalogueR::eQTLcatalogue_search_metadata(
            qtl_search = "Alasoo_2018",
            return_dataset_ids = FALSE,
            verbose = FALSE
        )
    )
    # IDs start with QTD
    expect_true(all(grepl("^QTD", result_ids)))
    # Labels contain Alasoo_2018
    expect_true(all(grepl("Alasoo_2018", result_labels)))
})
