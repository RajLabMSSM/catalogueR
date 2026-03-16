# Convert ENSEMBL IDs to HGNC gene symbols

Convert ENSEMBL IDs to HGNC gene symbols using the
[EnsDb.Hsapiens.v75](https://rdrr.io/pkg/EnsDb.Hsapiens.v75/man/package.html)
Bioconductor package.

## Usage

``` r
ensembl_to_hgnc(ensembl_ids, unique_only = TRUE, verbose = TRUE)
```

## Arguments

- ensembl_ids:

  Character vector of Ensembl gene IDs.

- unique_only:

  Only query unique entries `ensembl_ids`.

- verbose:

  Print messages.

## See also

Other utils:
[`hgnc_to_ensembl()`](https://rajlabmssm.github.io/catalogueR/reference/hgnc_to_ensembl.md),
[`merge_files()`](https://rajlabmssm.github.io/catalogueR/reference/merge_files.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ensembl_ids <- c("ENSG00000176697", "ENSG00000128573", "ENSG00000109743")
gene_symbols <- catalogueR::ensembl_to_hgnc(ensembl_ids = ensembl_ids)
} # }
```
