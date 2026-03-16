# Example colocalization results

Example colocalization results from running
[COLOC_run](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md).
on GWAS summary stats from all loci in Nalls23andMe_2019
([doi:10.1016/S1474-4422(19)30320-5](https://doi.org/10.1016/S1474-4422%2819%2930320-5)
). These results were published in [Schilder & Raj (Human Molecular
Genetics, 2021)](https://pubmed.ncbi.nlm.nih.gov/34617105/)

## Usage

``` r
COLOC_get_example_res(save_dir = tempdir(), full = FALSE)
```

## Source

` ##### Subset results #### # Pre-processing gwas.qtl_paths <- eQTLcatalogue_example_queries() coloc_QTLs <- COLOC_run(gwas.qtl_paths = gwas.qtl_paths, nThread = 4, top_snp_only = TRUE, save_path = "~/Desktop/coloc_results.tsv.gz") # Import pre-processed results URL <- file.path("https://github.com/RajLabMSSM/catalogueR/raw/master", "data/coloc_QTLs.rda") tmp <- file.path(tempdir(),basename(URL)) download.file(URL, tmp) piggyback::pb_upload(file = tmp, repo = "RajLabMSSM/catalogueR") ##### Full results #### URL <- file.path("https://github.com/RajLabMSSM/catalogueR/raw/master", "data/coloc_QTLs_full.rda") tmp <- file.path(tempdir(),basename(URL)) download.file(URL, tmp) piggyback::pb_upload(file = tmp, repo = "RajLabMSSM/catalogueR") `

## Arguments

- save_dir:

  Local directory to cache data in.

- full:

  Download the non-filtered version of the colocalization results
  (Default: `FALSE`).

## See also

Other coloc:
[`COLOC_corplot()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_corplot.md),
[`COLOC_get_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_res.md),
[`COLOC_heatmap()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_heatmap.md),
[`COLOC_merge_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_merge_res.md),
[`COLOC_report_summary()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_report_summary.md),
[`COLOC_run()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coloc_QTLs <- catalogueR::COLOC_get_example_res()
} # }
```
