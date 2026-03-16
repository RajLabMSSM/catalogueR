# List available eQTL datasets

Does some additional preprocessing of metadata to categorize tissue
types.

## Usage

``` r
eQTLcatalogue_list_datasets(
  save_dir = tempdir(),
  force_new = FALSE,
  include_imported = TRUE,
  verbose = FALSE
)
```

## Source

[eQTL-Catalogue GitHub
repo](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources)

## Arguments

- save_dir:

  Where to save the processed metadata.

- force_new:

  Re-import the metadata from GitHub instead of using metadata that
  comes included with catalogueR.

- include_imported:

  Include metadata from datasets that have not yet been fully
  re-processed by eQTL Catalogue's standardized pipeline.

- verbose:

  Print messages.

## Examples

``` r
if (FALSE) { # \dontrun{
meta <- catalogueR::eQTLcatalogue_list_datasets()
} # }
```
