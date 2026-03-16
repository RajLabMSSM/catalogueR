# Search eQTL Catalogue metadata

Searches through multiple relevant metadata columns to find eQTL
Catalogue datasets that match at least one of your substrings in a list.
All searches are case-insensitive. If `qtl_search=NULL`, will return all
available datasets.

## Usage

``` r
eQTLcatalogue_search_metadata(
  qtl_search = NULL,
  return_dataset_ids = TRUE,
  verbose = TRUE
)
```

## Arguments

- qtl_search:

  This function will automatically search for any datasets that match
  your criterion. For example, if you search "Alasoo_2018", it will
  query the datasets:

  - Alasoo_2018.macrophage_naive

  - Alasoo_2018.macrophage_Salmonella

  - Alasoo_2018.macrophage_IFNg+Salmonella

  You can be more specific about which datasets you want to include, for
  example by searching: "Alasoo_2018.macrophage_IFNg". You can even
  search by tissue or condition type (e.g. `c("blood","brain")`) and any
  QTL datasets containing those substrings (case-insensitive) in their
  name or metadata will be queried too.

- return_dataset_ids:

  Return dataset IDs (default: `TRUE`) instead of human-readable dataset
  labels.

- verbose:

  Show more (`=TRUE`) or fewer (`=FALSE`) messages.

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
qtl_datasets <- eQTLcatalogue_search_metadata(qtl_search = c(
    "Alasoo_2018",
    "monocyte"
))
qtl_datasets.brain <- eQTLcatalogue_search_metadata(qtl_search = "brain")
} # }
```
