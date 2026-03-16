# Prepare data for coloc plot

Use the coloc results (coloc_QTLs) to determine which full summary stats
(gwas.qtl_paths) to plot.

## Usage

``` r
COLOC_merge_res(
  gwas.qtl_paths,
  coloc_QTLs = NULL,
  qtl_thresh = 1e-05,
  coloc_thresh = 0.8,
  gwas_label = "GWAS",
  remove_extra_panes = TRUE,
  verbose = TRUE
)
```

## Arguments

- gwas.qtl_paths:

  Query results paths from
  [eQTLcatalogue_query](https://rajlabmssm.github.io/catalogueR/reference/eQTLcatalogue_query.md).

- coloc_QTLs:

  Colocalization results from
  [COLOC_run](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md).

- qtl_thresh:

  QTL uncorrected p-value ("pvalues.QTL") threshold.

- coloc_thresh:

  Colocalization Posterior Probability threshold, using the formula:
  `(PP.H3 + PP.H4 >= coloc_thresh) & (PP.H4 / PP.H3 >= 2)`.

- gwas_label:

  Label for the GWAS subplot.

- remove_extra_panes:

  Remove SNPs from non-significant panes.

- verbose:

  Print messages.

## See also

Other coloc:
[`COLOC_corplot()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_corplot.md),
[`COLOC_get_example_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_example_res.md),
[`COLOC_get_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_res.md),
[`COLOC_heatmap()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_heatmap.md),
[`COLOC_report_summary()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_report_summary.md),
[`COLOC_run()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md)
