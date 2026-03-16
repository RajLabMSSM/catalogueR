# eQTL Catalogue dataset metadata

List of all queryable tabix-indexed eQTL Catalogue datasets and their
associated systems/tissues/cell types.

## Usage

``` r
meta
```

## Format

An object of class `data.table` (inherits from `data.frame`) with 758
rows and 14 columns.

## Source

` meta <- eQTLcatalogue_list_datasets(force_new = TRUE) usethis::use_data(meta, overwrite = TRUE) `

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md)
