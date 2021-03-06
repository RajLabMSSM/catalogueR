% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{meta}
\alias{meta}
\title{eQTL Catalogue dataset metadata}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 397 rows and 12 columns.
}
\usage{
meta
}
\description{
List of all queryable tabix-indexed eQTL Catalogue datasets
and their associated systems/tissues/cell types.
}
\examples{
\dontrun{
meta <- eQTL_Catalogue.list_datasets(force_new=T)
# Some paths in the metadata were originally wrong. Has since been corrected by authors.
### meta <- meta \%>\% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path)) 
usethis::use_data(meta, overwrite=T)
}
}
\seealso{
Other eQTL Catalogue: 
\code{\link{BST1__Alasoo_2018.macrophage_IFNg}},
\code{\link{LRRK2__Alasoo_2018.macrophage_IFNg}},
\code{\link{MEX3C__Alasoo_2018.macrophage_IFNg}},
\code{\link{eQTL_Catalogue.fetch}()},
\code{\link{eQTL_Catalogue.header}},
\code{\link{eQTL_Catalogue.iterate_fetch}()},
\code{\link{eQTL_Catalogue.query}()},
\code{\link{eQTL_Catalogue.search_metadata}()},
\code{\link{fetch_restAPI}()},
\code{\link{fetch_tabix}()},
\code{\link{merge_gwas_qtl}()}
}
\concept{eQTL Catalogue}
\keyword{datasets}
