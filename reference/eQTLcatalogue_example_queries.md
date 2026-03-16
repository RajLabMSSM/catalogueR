# Paths to example eQTL Catalogue query results

Returns the paths to eQTL Catalogue query results stored within
catalogueR. Each file is a merged
[data.table](https://rdrr.io/pkg/data.table/man/data.table.html) of the
GWAS summary stats used to make the query, and the eQTL Catalogue query
results (which can contain data for multiple eGenes).

## Usage

``` r
eQTLcatalogue_example_queries(
  save_dir = tempdir(),
  fnames = c("BST1__Alasoo_2018.macrophage_IFNg+Salmonella.tsv.gz",
    "BST1__Alasoo_2018.macrophage_IFNg.tsv.gz",
    "BST1__Alasoo_2018.macrophage_naive.tsv.gz",
    "BST1__Alasoo_2018.macrophage_Salmonella.tsv.gz",
    "LRRK2__Alasoo_2018.macrophage_IFNg.tsv.gz",
    "LRRK2__Alasoo_2018.macrophage_naive.tsv.gz",
    "MEX3C__Alasoo_2018.macrophage_IFNg.tsv.gz",
    "MEX3C__Alasoo_2018.macrophage_naive.tsv.gz")
)
```

## Source

[doi:10.1016/S1474-4422(19)30320-5](https://doi.org/10.1016/S1474-4422%2819%2930320-5)

` paths <- catalogueR:::eQTLcatalogue_example_queries_local() fnames <- lapply(paths, function(x){ piggyback::pb_upload(file = x, repo = "RajLabMSSM/catalogueR") return(basename(x)) }) `

## Arguments

- save_dir:

  Local directory to cache data in.

- fnames:

  Character vector of file names to download.

## Value

Path to merged example GWAS-QTL summary statistics.

## Details

GWAS data originally comes from the Parkinson's disease GWAS by Nalls et
al., 2019 (The Lancet Neurology)
([doi:10.1016/S1474-4422(19)30320-5](https://doi.org/10.1016/S1474-4422%2819%2930320-5)
).

## Examples

``` r
if (FALSE) { # \dontrun{
gwas.qtl_paths <- catalogueR::eQTLcatalogue_example_queries()
} # }
```
