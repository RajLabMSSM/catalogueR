# Annotate QTL datasets with metadata

Annotate QTL datasets from *eQTL Catalogue* with metadata.

## Usage

``` r
eQTLcatalogue_annotate_tissues(dat, add_tissue_counts = FALSE)
```

## Arguments

- dat:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) of
  query results from
  [eQTLcatalogue_example_queries](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_example_queries.md).

- add_tissue_counts:

  Add a new column "Tissue_count" summarizing the number of datasets per
  tissue.

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- catalogueR::eQTLcatalogue_example_queries(
    fnames = "BST1__Alasoo_2018.macrophage_IFNg+Salmonella.tsv.gz")
dat <- data.table::fread(paths[1], nThread = 1)
dat_annot <- catalogueR::eQTLcatalogue_annotate_tissues(dat = dat)
} # }
```
