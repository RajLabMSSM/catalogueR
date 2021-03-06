% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coloc.R
\name{run_coloc}
\alias{run_coloc}
\title{Iteratively run coloc on merged GWAS-QTL datatables}
\usage{
run_coloc(
  gwas.qtl_paths,
  save_path = "./coloc_results.tsv.gz",
  nThread = 4,
  top_snp_only = T,
  split_by_group = F,
  method = "abf",
  PP_threshold = 0.8
)
}
\value{
If \code{top_snp_only=T}, returns SNP-level stats for only the SNP 
with the highest colocalization probability (\emph{SNP.PP.H4})
If \code{top_snp_only=T}, returns SNP-level stats for every SNP.
In either case, summary-level coloc stats are added in the columns
\emph{PP.H0}, \emph{PP.H1}, \emph{PP.H2}, \emph{PP.H3}, \emph{PP.H4}.
}
\description{
Runs colocalization tests (\code{\link{coloc::coloc.abf}}) on merged GWAS-QTL datatables
generated by \code{\link{catalogueR::eQTL_Catalogue.query}}.
}
\details{
Iterately runs coloc across each:
\itemize{
\item{QTL dataset}
\item{GWAS locus}
\item{QTL gene}
}
}
\examples{
# With built-in data
gwas.qtl_paths <- example_eQTL_Catalogue_query_paths() 
coloc_QTLs <- run_coloc(gwas.qtl_paths=gwas.qtl_paths, nThread=4, top_snp_only=T, save_path="~/Desktop/coloc_results.tsv.gz")


\dontrun{ 
# With full  Nalls et al data (not included)
gwas.qtl_paths <- list.files("/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019", recursive=T, full.names = T)
coloc_QTLs.Nalls2019 <- run_coloc(gwas.qtl_paths=gwas.qtl_paths[1:100], nThread=4, top_snp_only=T, save_path="~/Desktop/Nall2019.coloc_results.tsv.gz")
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
\code{\link{retrieve_sumstats_info}()}
}
\concept{coloc}
