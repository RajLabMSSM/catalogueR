# Merge files from a list of paths

Merge a list of files into one by stacking them on top of each other
(i.e. `rbind`).

## Usage

``` r
merge_files(file_paths, nThread = 1, verbose = TRUE)
```

## Arguments

- file_paths:

  Paths of files to import and merge into one
  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html).

- nThread:

  Number of threads to parallelize reading in files across.

- verbose:

  Print messages.

## See also

Other utils:
[`ensembl_to_hgnc()`](https://rajlabmssm.github.io/catalogueR/reference/ensembl_to_hgnc.md),
[`hgnc_to_ensembl()`](https://rajlabmssm.github.io/catalogueR/reference/hgnc_to_ensembl.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sumstats_paths <- echodata::get_Nalls2019_loci(limit_snps = 5)
merged_dat <- merge_files(file_paths = sumstats_paths)
} # }
```
