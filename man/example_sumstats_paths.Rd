% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{example_sumstats_paths}
\alias{example_sumstats_paths}
\title{Paths to example summary stats}
\source{
\url{https://www.biorxiv.org/content/10.1101/388165v3}
}
\usage{
example_sumstats_paths(Rlib_path = NULL)
}
\arguments{
\item{Rlib_path}{This function will automatically find your Rlib path, 
but you can override this by supplying it manually.}
}
\description{
Returns the paths to summary stats stored within \emph{catalogueR}.
Each file is the output of a locus that has been fine-mapping using \emph{echolocatoR}.
Data originally comes from the Parkinson's disease GWAS
by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
}
\details{
\describe{
  \item{SNP}{SNP RSID}
  \item{CHR}{Chromosome}
  \item{POS}{Genomic positiion (in basepairs)}
  ...
}
}
\examples{
sumstats_paths <- example_sumstats_paths()
}
\seealso{
Other Nalls23andMe_2019
These example files can be used: 
\code{\link{example_eQTL_Catalogue_query_paths}()}
}
\concept{Nalls23andMe_2019
These example files can be used}
