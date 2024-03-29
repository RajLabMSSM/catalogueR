% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COLOC_get_res.R
\name{COLOC_get_res}
\alias{COLOC_get_res}
\title{Run coloc on GWAS-QTL object}
\source{
\code{
library(dplyr)
paths <- catalogueR::eQTLcatalogue_example_queries()
query <- paths["BST1__Alasoo_2018.macrophage_IFNg"]

qtl.egene <- data.frame(
    query
    )[, grep("*.QTL$|qtl_id|SNP",
             colnames(query), value = TRUE)]
sorted_egenes <- qtl.egene |>
    dplyr::group_by(gene.QTL) |>
    dplyr::summarise(mean.P = mean(pvalue.QTL), min.P = min(pvalue.QTL)) |>
    dplyr::arrange(min.P)
qtl.egene <- subset(qtl.egene, gene.QTL == sorted_egenes$gene.QTL[1])

gwas.region <- data.frame(
    query
    )[, grep("*.QTL$|qtl_id",
             colnames(query),
             value = TRUE, invert = TRUE)]
#### Run ####
coloc_res <- catalogueR:::COLOC_get_res(qtl.egene = qtl.egene,
                                     gwas.region = gwas.region)
}
}
\usage{
COLOC_get_res(
  qtl.egene,
  gwas.region,
  merge_by_rsid = TRUE,
  coloc_thresh = 0.8,
  method = "abf",
  verbose = TRUE
)
}
\description{
Run coloc on GWAS-QTL object.
}
\seealso{
Other coloc: 
\code{\link{COLOC_corplot}()},
\code{\link{COLOC_get_example_res}()},
\code{\link{COLOC_heatmap}()},
\code{\link{COLOC_merge_res}()},
\code{\link{COLOC_report_summary}()},
\code{\link{COLOC_run}()}
}
\concept{coloc}
\keyword{internal}
