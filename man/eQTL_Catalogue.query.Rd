% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{eQTL_Catalogue.query}
\alias{eQTL_Catalogue.query}
\title{Iterate queries to \emph{eQTL Catalogue}}
\usage{
eQTL_Catalogue.query(
  sumstats_paths = NULL,
  output_dir = "./catalogueR_queries",
  qtl_search = NULL,
  use_tabix = T,
  conda_env = "echoR",
  nThread = 4,
  quant_method = "ge",
  infer_region = T,
  split_files = T,
  merge_with_gwas = T,
  force_new_subset = F,
  genome_build = "hg19",
  progress_bar = T,
  verbose = T
)
}
\arguments{
\item{sumstats_paths}{A list of paths to any number of summary stats files 
whose coordinates you want to use to make queries to eQTL Catalogue.  
If you wish to add custom names to the loci, simply add these as the names of the path list
 (e.g. \code{c(BST1="<path>/<to>/<BST1_file>", LRRK2="<path>/<to>/<LRRK2_file>")}). 
 Otherwise, loci will automatically named based on their min/max genomic coordinates. 
 
The minimum columns in these files required to make queries include: 
\describe{
\item{SNP}{RSID of each SNP.} 
\item{CHR}{Chromosome (can be in "chr12" or "12" format).}  
\item{POS}{Genomic position of each SNP.}  
\item{...}{Optional extra columns.}  
}}

\item{output_dir}{The folder you want the merged gwas/qtl results to be saved to 
(set \code{output_dir=F} if you don't want to save the results).  
If \code{split_files=F}, all query results will be merged into one and saved as \emph{<output_dir>/eQTL_Catalogue.tsv.gz}.  
If \code{split_files=T}, all query results will instead be split into smaller files and stored in \emph{<output_dir>/}.}

\item{qtl_search}{This function will automatically search for any datasets that match your criterion.
For example, if you search "Alasoo_2018", it will query the datasets:  
\itemize{
\item{Alasoo_2018.macrophage_naive}  
\item{Alasoo_2018.macrophage_Salmonella}  
\item{Alasoo_2018.macrophage_IFNg+Salmonella}
} 
You can be more specific about which datasets you want to include, 
for example by searching: "Alasoo_2018.macrophage_IFNg".
You can even search by tissue or condition type (e.g. \code{c("blood","brain")}) 
and any QTL datasets containing those substrings (case-insensitive) in their name or metadata will be queried too.}

\item{use_tabix}{Tabix is about ~17x faster (\emph{default:} =T) than the REST API (\emph{=F}).}

\item{nThread}{The number of CPU cores you want to use to speed up your queries through parallelization.}

\item{quant_method}{eQTL Catalogue actually contains more than just eQTL data.
For each dataset, the following kinds of QTLs can be queried:  
\describe{
\item{gene expression QTL}{\code{quant_method="ge"} (\emph{default}) or \code{quant_method="microarray"}, depending on the dataset. \strong{catalogueR} will automatically select whichever option is available.}  
\item{exon expression QTL}{\emph{*under construction*}  \code{quant_method="ex"}}
\item{transcript usage QTL}{\emph{*under construction*}  \code{quant_method="tx"}}
\item{promoter, splice junction and 3' end usage QTL}{\emph{*under construction*}  \code{quant_method="txrev"}}   
}}

\item{split_files}{Save the results as one file per QTL dataset (with all loci within each file).
If this is set to \code{=T}, then this function will return the list of paths where these files were saved.
A helper function is provided to import and merge them back together in R.  
If this is set to  \code{=F}, then this function will instead return one big merged data.table  
containing results from all QTL datasets and all loci.  
\code{=F} is not recommended when you have many large loci and/or many QTL datasets,  
because you can only fit so much data into memory.}

\item{merge_with_gwas}{Whether you want to merge your QTL query results with your GWAS data 
(convenient, but takes up more storage).}

\item{force_new_subset}{By default, \strong{catalogueR} will use any pre-existing files that match your query. 
Set \code{force_new_subset=T} to override this and force a new query.}

\item{genome_build}{The genome build of your query coordinates (e.g. \code{gwas_data}). 
If your coordinates are in \emph{hg19}, \strong{catalogueR} will automatically lift them over 
to \emph{hg38} (as this is the build that eQTL Catalogue uses).}

\item{progress_bar}{\code{progress_bar=T} allows progress to be monitored even when multithreading enabled.
Requires R package \code{\link{pbmcapply}}.}

\item{verbose}{Show more (\code{=T}) or fewer (\code{=F}) messages.}
}
\description{
Determines which datasets to query using \code{qtl_search}.
Uses coordinates from stored summary stats files (e.g. GWAS) 
to determine which regions to query from \emph{eQTL Catalogue}. 
Each locus file can be stored separately, 
or merged together to form one large file with all query results.
}
\examples{
sumstats_paths <- example_sumstats_paths()

# Merged results
# GWAS.QTL <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=1, force_new_subset=T, merge_with_gwas=F, progress_bar=T, split_files=F)
# Merged results (parallel)
GWAS.QTL <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=4, force_new_subset=T, merge_with_gwas=F, progress_bar=T, split_files=F)

# Split results
# gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=1, force_new_subset=T, merge_with_gwas=F, progress_bar=T) 
# Split results (parallel)
gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=4, force_new_subset=T, merge_with_gwas=F, progress_bar=T)
GWAS.QTL <- gather_files(file_paths = gwas.qtl_paths)

# Nalls et al example
\dontrun{
sumstats_paths_Nalls <- list.files("Fine_Mapping/Data/GWAS/Nalls23andMe_2019","Multi-finemap_results.txt", recursive = T, full.names = T)
names(sumstats_paths_Nalls) <- basename(dirname(dirname(sumstats_paths_Nalls)))
gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths_Nalls, output_dir="catalogueR_queries/Nalls23andMe_2019", merge_with_gwas=T, nThread=1, force_new_subset=T)
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
\code{\link{eQTL_Catalogue.search_metadata}()},
\code{\link{fetch_restAPI}()},
\code{\link{fetch_tabix}()},
\code{\link{merge_gwas_qtl}()},
\code{\link{meta}}
}
\concept{eQTL Catalogue}
