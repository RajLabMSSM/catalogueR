
#' Rapid querying, colocalization, and plotting of summary stats from the eQTL Catalogue
#' 
#' The functions in \strong{catalogueR} are partly derived from the following 
#' \href{http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html}{eQTL Catalogue tutorial}.
#' Additional eQTL Catalogue Resources:  
#' \itemize{
#' \item{\href{https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources }{GitHub}}  
#' \item{\href{https://www.ebi.ac.uk/eqtl/api-docs/}{In-depth API documentation}}  
#' \item{\href{ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv}{FTP server}}
#' } 
#' 
#' \strong{Notes on parallelization}: 
#' There's multiple levels to parallelize on.
#' You can only chooose one level at a time:  
#' \describe{
#' \item{\code{multithread_qtl=T}}{Across QTL datasets}  
#' \item{\code{multithread_loci=T}}{Across loci}
#' \item{\code{multithread_tabix=T}}{Within tabix files}  
#' }
#'  
#' You can also get a speedup by using tabix instead of the rest API
#' Test: For 3 loci, and X QTL datasets:
#' \describe{
#' \item{\strong{RESTful API:}}{7.5 minutes}  
#' \item{\strong{Tabix:}}{27 seconds (*clear winner! ~17x speedup)}
#' That said, if you're only query a small number of specific SNPs 
#' (rather than a large range), the RESTful API can sometimes be faster.  
#' }  
"_PACKAGE"
