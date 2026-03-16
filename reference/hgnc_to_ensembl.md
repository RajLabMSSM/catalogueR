# Convert HGNC gene symbols to ENSEMBL IDs

Convert HGNC gene symbols to ENSEMBL IDs using the
[EnsDb.Hsapiens.v75](https://rdrr.io/pkg/EnsDb.Hsapiens.v75/man/package.html)
Bioconductor package.

## Usage

``` r
hgnc_to_ensembl(gene_symbols, unique_only = TRUE, verbose = TRUE)
```

## Arguments

- gene_symbols:

  Character vector of HGNC gene IDs.

- unique_only:

  Only query unique entries `gene_symbols`.

- verbose:

  Print messages.

## See also

Other utils:
[`ensembl_to_hgnc()`](https://rajlabmssm.github.io/catalogueR/reference/ensembl_to_hgnc.md),
[`merge_files()`](https://rajlabmssm.github.io/catalogueR/reference/merge_files.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gene_symbols <- c("BDNF", "FOXP2", "BST1")
ensembl_ids <- catalogueR::hgnc_to_ensembl(gene_symbols)
} # }
```
