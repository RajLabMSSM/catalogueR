% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{retrieve_sumstats_info}
\alias{retrieve_sumstats_info}
\title{Merge coloc results with SNP-wise GWAS/QTL data}
\usage{
retrieve_sumstats_info(output_dir = "./catalogueR_queries", coloc_QTLs)
}
\description{
Merge coloc results with SNP-wise GWAS/QTL data
}
\examples{
\dontrun{
data("coloc_QTLs")
query_paths <- example_eQTL_Catalogue_query_paths()
output_dir <- unique(dirname(dirname(query_paths))) 
} 
}
\seealso{
Other coloc: 
\code{\link{COLOC.report_summary}()},
\code{\link{coloc_QTLs_full}},
\code{\link{coloc_QTLs}},
\code{\link{gather_colocalized_data}()},
\code{\link{get_colocs}()},
\code{\link{plot_coloc_summary}()},
\code{\link{run_coloc}()}
}
\concept{coloc}
