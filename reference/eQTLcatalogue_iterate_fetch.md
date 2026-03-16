# Iterate queries to *eQTL Catalogue*

Uses coordinates from stored summary stats files (e.g. GWAS) to
determine which regions to query from *eQTL Catalogue*.

## Usage

``` r
eQTLcatalogue_iterate_fetch(
  sumstats_paths,
  output_dir = file.path(tempdir(), "catalogueR_queries"),
  dataset_id,
  method = c("REST", "tabix"),
  quant_method = "ge",
  multithread_loci = TRUE,
  multithread_tabix = FALSE,
  split_files = TRUE,
  merge_with_gwas = FALSE,
  force_new_subset = FALSE,
  query_genome = "hg19",
  conda_env = "echoR_mini",
  nThread = 1,
  verbose = TRUE
)
```

## Source

` sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5) dataset_id <- catalogueR::eQTLcatalogue_list_datasets()$unique_label[1] GWAS.QTL <- catalogueR:::eQTLcatalogue_iterate_fetch( sumstats_paths = sumstats_paths, dataset_id = dataset_id, nThread = 1, split_files = FALSE) `

## Arguments

- sumstats_paths:

  A list of paths to any number of summary stats files whose coordinates
  you want to use to make queries to eQTL Catalogue. If you wish to add
  custom names to the loci, simply add these as the names of the path
  list (e.g.
  `c(BST1="<path>/<to>/<BST1_file>", LRRK2="<path>/<to>/<LRRK2_file>")`).
  Otherwise, loci will automatically named based on their min/max
  genomic coordinates.

  The minimum columns in these files required to make queries include:

  SNP

  :   RSID of each SNP.

  CHR

  :   Chromosome (can be in "chr12" or "12" format).

  POS

  :   Genomic position of each SNP.

  ...

  :   Optional extra columns.

- output_dir:

  The folder you want the merged gwas/qtl results to be saved to (set to
  `NULL` to not save the results). If `split_files=FALSE`, all query
  results will be merged into one and saved as
  *\<output_dir\>/eQTLcatalogue_tsv.gz*. If `split_files=TRUE`, all
  query results will instead be split into smaller files and stored in
  *\<output_dir\>/*.

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

- split_files:

  Save the results as one file per QTL dataset (with all loci within
  each file). If this is set to `=TRUE`, then this function will return
  the list of paths where these files were saved. A helper function is
  provided to import and merge them back together in R. If this is set
  to `=FALSE`, then this function will instead return one big merged
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  containing results from all QTL datasets and all loci. `=FALSE` is not
  recommended when you have many large loci and/or many QTL datasets,
  because you can only fit so much data into memory.

- merge_with_gwas:

  Whether you want to merge your QTL query results with your GWAS data
  (convenient, but takes up more storage).

- force_new_subset:

  By default, **catalogueR** will use any pre-existing files that match
  your query. Set `force_new_subset=T` to override this and force a new
  query.

- query_genome:

  The genome build of your query coordinates (e.g. `query_dat`). If your
  coordinates are in *hg19*, **catalogueR** will automatically lift them
  over to *hg38* (as this is the build that eQTL Catalogue uses).

- conda_env:

  Conda environment to search for tabix executable in.

- nThread:

  The number of CPU cores you want to use to speed up your queries
  through parallelization.

- verbose:

  Show more (`=TRUE`) or fewer (`=FALSE`) messages.

## See also

Other eQTL Catalogue:
[`eQTLcatalogue_fetch()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_fetch.md),
[`eQTLcatalogue_header`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_header.md),
[`eQTLcatalogue_query()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md),
[`eQTLcatalogue_search_metadata()`](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_search_metadata.md),
[`fetch_restAPI()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_restAPI.md),
[`fetch_tabix()`](https://rajlabmssm.github.io/catalogueR/reference/fetch_tabix.md),
[`merge_gwas_qtl()`](https://rajlabmssm.github.io/catalogueR/reference/merge_gwas_qtl.md),
[`meta`](https://rajlabmssm.github.io/catalogueR/reference/meta.md)
