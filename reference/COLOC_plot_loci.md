# Faceted Manhattan plots of QTL datasets

Use the coloc results (`coloc_QTLs`) to determine which full summary
statistics (`gwas.qtl_paths`) to plot.

## Usage

``` r
COLOC_plot_loci(
  gwas.qtl_paths = NULL,
  coloc_QTLs = NULL,
  plot_dat = NULL,
  qtl_thresh = 1e-05,
  coloc_thresh = 0.8,
  gwas_label = "GWAS",
  remove_extra_panes = TRUE,
  y_facet_angle = 0,
  x_facet_angle = 270,
  show_plot = TRUE,
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

- plot_dat:

  \[Optional\] Pre-computed plot data from
  [COLOC_merge_res](https://rajlabmssm.github.io/catalogueR/reference/COLOC_merge_res.md).

- qtl_thresh:

  QTL uncorrected p-value ("pvalues.QTL") threshold.

- coloc_thresh:

  Colocalization Posterior Probability threshold, using the formula:
  `(PP.H3 + PP.H4 >= coloc_thresh) & (PP.H4 / PP.H3 >= 2)`.

- gwas_label:

  Label for the GWAS subplot.

- remove_extra_panes:

  Remove SNPs from non-significant panes.

- y_facet_angle:

  Angle of the y-axis facet labels.

- x_facet_angle:

  Angle of the x-axis facet labels.

- show_plot:

  Whether to print the plot.

- verbose:

  Print messages.

## Examples

``` r
if (FALSE) { # \dontrun{
coloc_QTLs <- catalogueR::COLOC_get_example_res()
gwas.qtl_paths <- catalogueR::eQTLcatalogue_example_queries()
gg_gwas.qtl <- catalogueR::COLOC_plot_loci(
    gwas.qtl_paths = gwas.qtl_paths,
    coloc_QTLs = coloc_QTLs,
    coloc_thresh = .5,
    qtl_thresh = .005,
    remove_extra_panes = FALSE)
} # }
```
