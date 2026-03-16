# eQTL Catalogue tabix header

The eQTL Catalogue tabix header (column names) is not always included in
the file.

## Usage

``` r
data("eQTLcatalogue_header")
```

## Format

An object of class `character` of length 19.

## Source

` eQTLcatalogue_header <- eQTLcatalogue_get_header(force_new_header = TRUE) usethis::use_data(eQTLcatalogue_header, overwrite = TRUE) `

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)
