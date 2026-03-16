# Plot summary of coloc tests: heatmap

Plots the output from
[COLOC_run](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md).

## Usage

``` r
COLOC_heatmap(
  coloc_QTLs,
  coloc_thresh = 0.8,
  label_snp_groups = FALSE,
  save_dir = tempdir(),
  show_plot = TRUE,
  verbose = TRUE
)
```

## Arguments

- coloc_QTLs:

  Colocalization results from
  [COLOC_run](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md).

- coloc_thresh:

  Colocalization Posterior Probability threshold, using the formula:
  `(PP.H3 + PP.H4 >= coloc_thresh) & (PP.H4 / PP.H3 >= 2)`.

- label_snp_groups:

  Whether to label SNP groups.

- save_dir:

  Directory to save plot in.

- show_plot:

  Whether to print the plot.

- verbose:

  Print messages.

## See also

Other coloc:
[`COLOC_corplot()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_corplot.md),
[`COLOC_get_example_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_example_res.md),
[`COLOC_get_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_res.md),
[`COLOC_merge_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_merge_res.md),
[`COLOC_report_summary()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_report_summary.md),
[`COLOC_run()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coloc_QTLs <- catalogueR::COLOC_get_example_res()
gg_coloc <- catalogueR::COLOC_heatmap(
    coloc_QTLs = coloc_QTLs,
    coloc_thresh = .5)
} # }
```
