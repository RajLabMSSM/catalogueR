# 2. Query eQTL Catalogue datasets by region

2.2 Method 2: RESTful API Slower than tabix (unless you're only querying
several specific SNPs).

## Usage

``` r
fetch_restAPI(
  dataset_id,
  query_granges,
  quant_method = "ge",
  is_gwas = FALSE,
  size = NULL,
  timeout = 5 * 60,
  verbose = TRUE
)
```

## Source

` data("meta") query_granges <- echotabix::construct_query(query_dat = echodata::BST1) qtl.subset <- fetch_restAPI(dataset_id = meta$dataset_id[1], query_granges = query_granges) `

## Arguments

- quant_method:

  eQTL Catalogue actually contains more than just eQTL data. For each
  dataset, the following kinds of QTLs can be queried:

  gene expression QTL

  :   `quant_method="ge"` (*default*) or `quant_method="microarray"`,
      depending on the dataset. **catalogueR** will automatically select
      whichever option is available.

  exon expression QTL

  :   *\*under construction\** `quant_method="ex"`

  transcript usage QTL

  :   *\*under construction\** `quant_method="tx"`

  promoter, splice junction and 3' end usage QTL

  :   *\*under construction\** `quant_method="txrev"`

- verbose:

  Show more (`=TRUE`) or fewer (`=FALSE`) messages.

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)
