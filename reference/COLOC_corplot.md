# Plot summary of coloc tests: heatmap

Plots the output from
[COLOC_run](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md).

## Usage

``` r
COLOC_corplot(
  coloc_QTLs,
  coloc_thresh = 0.8,
  label_top_snps = 1,
  facets = Study + System + Tissue ~ Locus.eGene,
  show_plot = TRUE,
  save_dir = tempdir(),
  seed = 2022,
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

- label_top_snps:

  Label the top n SNPs per `facets` group, sorted by `SNP.PP.H4` from
  coloc.

- facets:

  Formula specifying the faceting variables for the plot. Passed to
  [facet_grid](https://ggplot2.tidyverse.org/reference/facet_grid.html).

- show_plot:

  Whether to print the plot.

- save_dir:

  Directory to save plot in.

- seed:

  Random seed for reproducibility of label placement via
  [geom_label_repel](https://ggrepel.slowkow.com/reference/geom_text_repel.html).

- verbose:

  Print messages.

## See also

Other coloc:
[`COLOC_get_example_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_example_res.md),
[`COLOC_get_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_res.md),
[`COLOC_heatmap()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_heatmap.md),
[`COLOC_merge_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_merge_res.md),
[`COLOC_report_summary()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_report_summary.md),
[`COLOC_run()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md)

## Examples

``` r
if (FALSE) { # \dontrun{
coloc_QTLs <- catalogueR::COLOC_get_example_res()
gg_coloc <- COLOC_corplot(coloc_QTLs = coloc_QTLs,
                          coloc_thresh = .5)
} # }
```
