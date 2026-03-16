# Run coloc on GWAS-QTL object

Run coloc on GWAS-QTL object.

## Usage

``` r
COLOC_get_res(
  qtl.egene,
  gwas.region,
  merge_by_rsid = TRUE,
  coloc_thresh = 0.8,
  method = "abf",
  verbose = TRUE
)
```

## Source

` library(dplyr) paths <- catalogueR::eQTLcatalogue_example_queries() query <- paths["BST1__Alasoo_2018.macrophage_IFNg"] qtl.egene <- data.frame( query )[, grep("*.QTL$|dataset_id|SNP", colnames(query), value = TRUE)] sorted_egenes <- qtl.egene |> dplyr::group_by(gene.QTL) |> dplyr::summarise(mean.P = mean(pvalue.QTL), min.P = min(pvalue.QTL)) |> dplyr::arrange(min.P) qtl.egene <- subset(qtl.egene, gene.QTL == sorted_egenes$gene.QTL[1]) gwas.region <- data.frame( query )[, grep("*.QTL$|dataset_id", colnames(query), value = TRUE, invert = TRUE)] #### Run #### coloc_res <- catalogueR:::COLOC_get_res(qtl.egene = qtl.egene, gwas.region = gwas.region) `

## See also

Other coloc:
[`COLOC_corplot()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_corplot.md),
[`COLOC_get_example_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_get_example_res.md),
[`COLOC_heatmap()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_heatmap.md),
[`COLOC_merge_res()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_merge_res.md),
[`COLOC_report_summary()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_report_summary.md),
[`COLOC_run()`](https://rajlabmssm.github.io/catalogueR/reference/COLOC_run.md)
