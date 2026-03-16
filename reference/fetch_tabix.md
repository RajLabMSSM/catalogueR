# Query eQTL Catalogue:tabix/echotabix

Query eQTL Catalogue datasets by region with tabix or echotabix. Faster
alternative to REST API.

## Usage

``` r
fetch_tabix(
  dataset_id,
  query_granges,
  query_method = c("rsamtools", "variantannotation", "conda", "seqminer"),
  quant_method = "ge",
  nThread = 1,
  conda_env = "echoR_mini",
  verbose = TRUE
)
```

## Source

[eQTL Catalogue blocking tabix
requests](https://github.com/RajLabMSSM/catalogueR/issues/5)

` data("meta"); query_granges <- echotabix::construct_query(query_dat = echodata::BST1) qtl.subset <- catalogueR:::fetch_tabix(unique_label=meta$unique_label[2], query_granges=query_granges) `

## Arguments

- query_granges:

  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object to be used for querying the `target_path` file. Alternatively,
  can be variant-level summary statistics to be converted into a
  [GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object by
  [construct_query](https://rdrr.io/pkg/echotabix/man/construct_query.html).

- query_method:

  Method used for querying. See
  [query](https://rdrr.io/pkg/echotabix/man/query.html) for available
  options.

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

- nThread:

  The number of CPU cores you want to use to speed up your queries
  through parallelization.

- conda_env:

  Conda environment to search for tabix executable in.

- verbose:

  Show more (`=TRUE`) or fewer (`=FALSE`) messages.

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)
