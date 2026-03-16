# Query eQTL Catalogue

Query eQTL Catalogue datasets with multiple methods options.

## Usage

``` r
eQTLcatalogue_fetch(
  dataset_id,
  query_granges,
  method = c("REST", "tabix"),
  quant_method = "ge",
  multithread_tabix = FALSE,
  add_dataset_id = TRUE,
  convert_genes = TRUE,
  suffix = ".QTL",
  timeout = 5 * 60,
  conda_env = "echoR_mini",
  nThread = 1,
  verbose = TRUE
)
```

## Arguments

- dataset_id:

  Unique eQTL Catalogue ID assigned in metadata ("dataset_id" column in
  `data(meta)`).

- query_granges:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) with
  at least columns `CHR` and `POS` defining the genomic region to query.

- method:

  Method for querying eQTL Catalogue:

  REST (default)

  :   Uses the REST API. Slow but can be used by anyone.

  tabix

  :   Uses tabix [query](https://rdrr.io/pkg/echotabix/man/query.html).
      Fast, but requires the user to first get their IP address
      whitelisted by the EMBL-EBI server admin by putting in a request
      [here](https://www.ebi.ac.uk/about/contact/support/).

  *Note*: "tabix" is about ~17x faster than the REST API, but is
  currently a far less reliable method than the REST API because tabix
  tends to get blocked by eQTL Catalogue's firewall. See
  [here](https://github.com/RajLabMSSM/catalogueR/issues/5) for more
  details.

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

- multithread_tabix:

  Multi-thread across within a single tabix file query (good when you
  have one-several large loci).

- add_dataset_id:

  Add "dataset_id" (i.e. "dataset_id") column to the query result.

- convert_genes:

  Convert Ensembl IDs to HGNC symbols.

- suffix:

  Character suffix to append to QTL column names (default: `".QTL"`).
  Set to `NULL` to skip.

- timeout:

  Timeout in seconds for the REST API query (default: `5*60`).

- conda_env:

  Conda environment to search for tabix executable in.

- nThread:

  The number of CPU cores you want to use to speed up your queries
  through parallelization.

- verbose:

  Show more (`=TRUE`) or fewer (`=FALSE`) messages.

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_iterate_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_iterate_fetch.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("meta")
query_granges <- echodata::BST1
GWAS.QTL_manual <- catalogueR::eQTLcatalogue_fetch(
  query_granges = query_granges,
  dataset_id = meta$dataset_id[1])
} # }
```
