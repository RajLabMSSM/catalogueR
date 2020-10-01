

# Get the header (according to email with Kaur Alassoo)
tabix_header <- function(tabix_path=NULL, 
                         force_new_header=F){
  if(is.null(tabix_path)){tabix_path <- "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_naive.all.tsv.gz"}
  header_path <- system.file("eQTL_Catalogue_resources","tabix_header.txt",package = "catalogueR")
  
  if(file.exists(header_path) & force_new_header==F){
    header <- as.list(data.table::fread(header_path, nThread = 1))[[1]]
  } else { 
    # Read in header data
    header <-  colnames(data.table::fread(cmd = paste("curl -s",tabix_path,"| zcat | head -n 1"), nThread = 1))
    # Save header data 
    dir.create(dirname(header_path), showWarnings = F, recursive = T)
    data.table::fwrite(list(header), file=header_path, nThread = 1)
  }
  return(header)
}



 
#' 2. Query eQTL Catalogue datasets by region
#' 
#' 2.1 Method 1: Tabix
#' Faster alternative to REST API.
#' @inheritParams eQTL_Catalogue.query
#' @family eQTL Catalogue 
#' @examples 
#' data("meta"); data("BST1");
#' qtl.subset <- fetch_tabix(unique_id=meta$unique_id[2], gwas_dat=BST1)
fetch_tabix <- function(unique_id,
                        quant_method="ge",
                        infer_region=T,
                        gwas_data=NULL,
                        chrom=NULL,
                        bp_lower=NULL,
                        bp_upper=NULL,
                        is_gwas=F,
                        nThread=4,
                        verbose=T){
  # quant_method="ge"; infer_region=T;is_gwas=F; remove_tmp=F;  add_chr=T 
  check_coord_input(gwas_data=gwas_data, 
                    chrom=chrom, 
                    bp_lower=bp_lower, bp_upper=bp_upper)
  tabix.start = Sys.time()
  # Get region
  if(infer_region & !is.null(gwas_data)){
    print("++ Inferring coordinates from gwas_data", v=verbose)
    chrom <- unique(gwas_data$CHR)
    if(length(chrom)>1){stop("More than one chromosome detected.")}
    bp_lower <- min(gwas_data$POS)
    bp_upper <- max(gwas_data$POS)
  }
  region <- paste0(chrom,":",bp_lower,"-",bp_upper)
  # printer("+ TABIX:: Querying region:", region) 
  meta.sub <- choose_quant_method(ui=unique_id, 
                                  qm=quant_method, 
                                  verbose=verbose) 
  header <- tryCatch(expr = {
    tabix_header(tabix_path = meta.sub$ftp_path, 
                 force_new_header = F)
  },error= function(e){
    tabix_header(tabix_path = meta.sub$ftp_path,
                 force_new_header = F)
 })
  
  # Run tabix
  # Read directly into R rather than saving tabix subset
  qtl.subset <- data.table::fread(cmd=paste("tabix",
                                            # "--print-header",
                                            meta.sub$ftp_path,
                                            region),
                                  nThread = nThread)
  if(length(header) != ncol(qtl.subset)){
    header <- tabix_header(tabix_path = meta.sub$ftp_path, 
                           force_new_header = T)
  }
  colnames(qtl.subset) <- paste0(header,".QTL")
  tabix.end = Sys.time() 
  printer("eQTL_Catalogue::",nrow(qtl.subset),"SNPs returned in", round(as.numeric(tabix.end-tabix.start),1),"seconds.", v=verbose)  
  return(qtl.subset)
}


 

# Sub-function of fetch_restAPI
fetch_from_eqtl_cat_API <- function(link,
                                    is_gwas = F){
  nullToNA <- function(x) {
    x[sapply(x, is.null)] <- NA
    return(x)
  }
  if (is_gwas) {
    cols_to_nest <- c("variant_id", "chromosome", "base_pair_location", "trait",
                      "p_value", "ci_lower", "ci_upper", "beta",
                      "effect_allele", "other_allele", "effect_allele_frequency",
                      "odds_ratio", "study_accession", "code")
  } else {
    cols_to_nest <- c("study_id", "qtl_group", "rsid",
                      "chromosome", "position", "pvalue", "condition_label",
                      "tissue_label", "molecular_trait_id", "gene_id", "ac",
                      "ref", "beta",  "variant", "an", "median_tpm",  "condition",
                      "r2", "alt", "type", "maf",  "tissue")
  }
  is_paginated <- !stringr::str_detect(link,"paginate=False")
  # message("isPagined:", is_paginated)
  page = 1
  merged_summaries <- data.frame()
  while(!is.null(link)){
    # print(paste0("Fetching page #",page))
    api_raw_data <- jsonlite::fromJSON(link, simplifyDataFrame = T, flatten = T)
    link <- api_raw_data$`_links`$`next`$href
    if (rlang::is_empty(api_raw_data$`_embedded`$associations)) {
      return(merged_summaries)
    }
    eqtl_raw_list_data <- do.call(rbind, lapply(api_raw_data$`_embedded`$associations, rbind))
    eqtl_data <- nullToNA(eqtl_raw_list_data) %>% as.matrix() %>% dplyr::as_tibble()
    if (is_paginated) { eqtl_data <- dplyr::select(eqtl_data, -c("_links")) }
    eqtl_data <- tidyr::unnest(eqtl_data, cols = cols_to_nest)
    if (!is.null(link)) {
      page <- page + 1
    }
    merged_summaries <- merged_summaries %>% rbind(eqtl_data)
  }
  return(merged_summaries[cols_to_nest])
}




#' 2. Query eQTL Catalogue datasets by region
#' 
#' 2.2 Method 2: RESTful API
#' Slower than tabix (unless you're only querying several specific SNPs).
#' @inheritParams eQTL_Catalogue.query
#' @family eQTL Catalogue
#' @examples 
#' data("meta"); data("BST1");
#' qtl.subset <- fetch_restAPI(unique_id=meta$unique_id[1], gwas_data=BST1)
fetch_restAPI <- function(unique_id, #Alasoo_2018.macrophage_naive
                          quant_method="ge",
                          infer_region=T,
                          gwas_data=NULL,
                          chrom=NULL,
                          bp_lower=NULL,
                          bp_upper=NULL,
                          is_gwas=F, # refers to the datasets being queried
                          size=NULL,
                          verbose=T) {
  restAPI.start = Sys.time()
  # Get region
  if(infer_region & !is.null(gwas_data)){
    print("+ Inferring coordinates from gwas_data", v=verbose)
    chrom <- gsub("chr","",unique(gwas_data$CHR)[1])
    bp_lower <- min(gwas_data$POS)
    bp_upper <- max(gwas_data$POS)
  } 
  meta.sub <- choose_quant_method(ui=unique_id, 
                                  qm=quant_method, 
                                  verbose=verbose) 
  # gene_id <- hgnc_to_ensembl(gene_symbols = "LRRK2")
  link <- paste0("http://www.ebi.ac.uk/eqtl/api/",
                 "chromosomes/",chrom,
                 "/associations?paginate=False",
                 # Study name
                 "&study=",meta.sub$study,
                 # Condition
                 "&qtl_group=",meta.sub$qtl_group,
                 # ENSEMBL gene id
                 # "&gene_id=",gene_id,
                 # ENSEMBL molecular trait id
                 # "&molecular_trait_id=",molecular_trait_id,
                 # gene expression, transcript, etc.
                 "&quant_method=",quant_method,
                 # genomic position limits
                 "&bp_lower=",bp_lower,
                 "&bp_upper=",bp_upper,
                 if(!is.null(size)) paste0("&size=",size) else "")
  
  message("+ eQTL_Catalogue:: Fetching: ")
  message(link)
  qtl.subset <- fetch_from_eqtl_cat_API(link = link,  is_gwas = is_gwas)
  colnames(qtl.subset) <-  paste0(colnames(qtl.subset),".QTL")  
  restAPI.end = Sys.time()
  printer("+ eQTL_Catalogue::",nrow(qtl.subset),"SNPs returned in", round(as.numeric(restAPI.end-restAPI.start),1),"seconds.", v=verbose)  
  return(data.table::data.table(qtl.subset))
}





#' 2. Query eQTL Catalogue datasets by region
#'
#' Choose between tabix (faster for larger queries) 
#' or RESTful API (faster for small queries).
#' @inheritParams eQTL_Catalogue.query
#' @param bp_lower Minimum basepair position of the query window.
#' @param bp_upper Maxmimum basepair position of the query window.
#' @family eQTL Catalogue
#' @examples 
#' data("meta"); data("BST1");
#' gwas.qtl <- eQTL_Catalogue.fetch(unique_id=meta$unique_id[1], gwas_data=BST1)
eQTL_Catalogue.fetch <- function(unique_id,
                                 quant_method="ge",
                                 infer_region=T,
                                 gwas_data=NULL,
                                 is_gwas=F,
                                 nThread=1,
                                 use_tabix=T,
                                 chrom=NULL,
                                 bp_lower=NULL,
                                 bp_upper=NULL, 
                                 multithread_tabix=F,
                                 add_qtl_id=T,
                                 convert_genes=T,
                                 verbose=T){
  if(use_tabix){
    # Tabix is about ~17x faster than the REST API.
    gwas.qtl <- fetch_tabix(unique_id=unique_id,
                            quant_method=quant_method,
                            infer_region=T,
                            gwas_data=gwas_data,
                            chrom=chrom, 
                            bp_lower=bp_lower, 
                            bp_upper=bp_upper,
                            is_gwas=F,
                            nThread=if(multithread_tabix) nThread else 1, 
                            verbose=verbose)
  } else {
    gwas.qtl <- fetch_restAPI(unique_id=unique_id,
                              quant_method=quant_method,
                              infer_region=T,
                              gwas_data=gwas_data,
                              chrom=chrom, 
                              bp_lower=bp_lower, 
                              bp_upper=bp_upper,
                              is_gwas=F, 
                              verbose=verbose)
  }
  # Post=processing
  if(add_qtl_id){
    gwas.qtl <- data.table::data.table(qtl_id=unique_id, gwas.qtl)
  } 
  # Convert genes
  if(convert_genes){
    try({
      gene_dict <- ensembl_to_hgnc(ensembl_ids = gwas.qtl$gene_id.QTL,  
                                   verbose = verbose)
      gwas.qtl$gene.QTL <- gene_dict[gwas.qtl$molecular_trait_object_id.QTL]
    }) 
  } 
  return(gwas.qtl)
}





#' Merge GWAS data (query) and QTL data (results)
#' 
#' @family eQTL Catalogue
#' @inheritParams eQTL_Catalogue.query
merge_gwas_qtl <- function(gwas_data,
                           qtl.subset, 
                           verbose=T){ 
  printer("++ Merging GWAS data and QTL query results.",v=verbose)
  gwas.qtl <- tryCatch(expr = {
    # Merging and allele flipping
    gwas.qtl <- data.table::merge.data.table(x = data.table::data.table(gwas_data),
                                             y = data.table::data.table(qtl.subset),
                                             # all.x = T,
                                             by.x = c("SNP"), # effect_allele
                                             by.y = c("rsid.QTL") ) %>%
      # subset(effect.is.ref|effect.is.alt) %>%
      data.table::data.table()
    if("A1" %in% colnames(gwas.qtl) & "A2" %in% colnames(gwas.qtl)){
      gwas.qtl <- gwas.qtl %>% dplyr::mutate(effect.is.ref=ifelse(A1==ref.QTL,T,F),
                                             effect.is.alt=ifelse(A2==alt.QTL,T,F) )
    }  
    return(gwas.qtl)
  }, error=function(e){return(qtl.subset)})
 return(gwas.qtl)
}



 
#' Iterate queries to \emph{eQTL Catalogue}
#' 
#' Uses coordinates from stored summary stats files (e.g. GWAS) 
#' to determine which regions to query from \emph{eQTL Catalogue}. 
#' @inheritParams eQTL_Catalogue.query
#' @param progress_bar Show progress bar during parallelization across loci.
#'  \emph{WARNING!}: Progress bar (via \code{\link{pbmclapply}}) only works
#'   on Linux/Unix systems (e.g. mac) and NOT on Windows.
#' @family eQTL Catalogue 
#' @examples 
#' sumstats_paths <- example_sumstats_paths()
#' qtl_id  <- eQTL_Catalogue.list_datasets()$unique_id[1]
#' GWAS.QTL <- eQTL_Catalogue.iterate_fetch(sumstats_paths=sumstats_paths, qtl_id=qtl_id, force_new_subset=T, multithread_loci=T, nThread=1, split_files=F, progress_bar=F)
eQTL_Catalogue.iterate_fetch <- function(sumstats_paths,  
                                         output_dir="./catalogueR_queries",
                                         qtl_id,
                                         quant_method="ge",
                                         infer_region=T, 
                                         use_tabix=T,
                                         multithread_loci=T,
                                         multithread_tabix=F,
                                         nThread=4,
                                         split_files=T,
                                         merge_with_gwas=F,
                                         force_new_subset=F,
                                         progress_bar=F,
                                         genome_build="hg19",
                                         verbose=T){
  # WARNING!: pbmclapply only worked on Linux/Unix systems (e.g. mac) and NOT on Windows.
  lapply_func <- progress_bar_check(progress_bar=progress_bar,
                                    verbose=verbose)
  ##########-----# ITERATE ACROSS LOCI ---------------------
  GWAS.QTL <- lapply_func(1:length(sumstats_paths), function(i,
                                                             .sumstats_paths=sumstats_paths,
                                                             .qtl_id=qtl_id,
                                                             .quant_method=quant_method,
                                                             .infer_region=infer_region,  
                                                             .nThread=nThread,
                                                             .multithread_tabix=multithread_tabix,
                                                             .use_tabix=use_tabix,
                                                             .force_new_subset=force_new_subset,
                                                             .genome_build=genome_build,
                                                             .split_files=split_files,
                                                             .merge_with_gwas=merge_with_gwas,
                                                             .verbose=verbose){
    # Import GWAS data
    # Have to do it this way in order to get names of sumstats_paths
    loc_path <- .sumstats_paths[i]
    gwas_data <- data.table::fread(loc_path, nThread = 1)
    gwas_data$CHR <- gsub("chr","",gwas_data$CHR)   # get rid of "chr" just in case
    
    # Name Locus
    if(!is.null(loc_path)){
      printer("++ Extracting locus name from `sumstats_paths` names.", v=.verbose)
      loc <- names(loc_path)
      if(!"Locus" %in% colnames(gwas_data)){gwas_data <- cbind(Locus=loc, gwas_data)}
    } else {
      if("Locus" %in% colnames(gwas_data)){
        printer("++ Extracting locus name from GWAS file.", v=.verbose)
        loc <- unique(gwas_data$Locus)[1]
      } else { 
        loc <- construct_locus_name(gwas_data, verbose=.verbose) 
        gwas_data <- cbind(Locus=loc, gwas_data)
      }
    } 
    
    message("_+_+_+_+_+_+_+_+_--- Locus: ",loc," ---_+_+_+_+_+_+_+_+_")
    # Test if query file already exists
    split_path <- make_split_path(output_dir=output_dir, 
                                  qtl_id=.qtl_id, 
                                  loc=loc) 
    if(file.exists(split_path) & .force_new_subset==F){
      printer("++ Using pre-existing file...", v=.verbose)
      qtl.subset <- data.table::fread(split_path, nThread = 1)
    } else {
      
      # Convert from GRCh37 to GRCh38
      if(.genome_build %in% c('hg19','hg18')){
        gr.lifted <- liftover(gwas_data = gwas_data, 
                              build.conversion = paste0(.genome_build,".to.hg38"), 
                              verbose=F)
        gwas_data <- data.frame(gr.lifted)  %>% 
          dplyr::mutate(CHR=gsub("chr","",seqnames),
                        POS=start) %>% 
          dplyr::select(-c("seqnames","start","end","width","strand")) %>%
          data.table::as.data.table()
      }
       
      qtl.subset <- data.table::data.table()
      qtl.subset <- tryCatch(expr = {
        eQTL_Catalogue.fetch(unique_id=.qtl_id,
                             quant_method=.quant_method,
                             infer_region=.infer_region,
                             gwas_data=gwas_data,
                             is_gwas=F,
                             nThread=.nThread,
                             multithread_tabix=.multithread_tabix,
                             use_tabix=.use_tabix,
                             chrom=NULL,
                             bp_upper=NULL,
                             bp_lower=NULL,
                             add_qtl_id=T,
                             convert_genes=T,
                             verbose = .verbose)
        },
        error=function(x){data.table::data.table()})
      # Merge results
      if(.merge_with_gwas){
        gwas.qtl <- tryCatch(expr = {
          merge_gwas_qtl(gwas_data=gwas_data, 
                         qtl.subset=qtl.subset, 
                         verbose=.verbose)
        }, 
        error=function(e){data.table::data.table()}) 
      } else {
        gwas.qtl <- qtl.subset
      }  
      tryCatch({
        # Add locus name
        printer("++ Adding `Locus.GWAS` column.", v=.verbose)
        gwas.qtl <- cbind(Locus.GWAS=loc,  
                          gwas.qtl)
        # Save
        if(.split_files){
          printer("++ Saving split file ==>",split_path, v=.verbose)
          dir.create(dirname(split_path), showWarnings = F, recursive = T)
          data.table::fwrite(gwas.qtl, split_path, sep="\t", nThread = 1)
        }
        # Return
        if(.split_files){return(split_path)} else {return(gwas.qtl)} 
      }, error=function(e){return(NULL)})
      
    }  
    
  }, mc.cores = if(multithread_loci) nThread else 1) ## END ITERATE ACROSS LOCI
  
  # Return
  if(split_files){
    printer("+ Returning list of split query results files.", v=verbose)
    return(unlist(GWAS.QTL))
    } else {
    printer("+ Returning merged data.table of query results.", v=verbose)
    GWAS.QTL <- data.table::rbindlist(GWAS.QTL, fill = T)  
    return(GWAS.QTL)
  } 
  message(" ")
}





#' Iterate queries to \emph{eQTL Catalogue}
#' 
#' Determines which datasets to query using \code{qtl_search}.
#' Uses coordinates from stored summary stats files (e.g. GWAS) 
#' to determine which regions to query from \emph{eQTL Catalogue}. 
#' Each locus file can be stored separately, 
#' or merged together to form one large file with all query results.
#' 
#' @param sumstats_paths A list of paths to any number of summary stats files 
#' whose coordinates you want to use to make queries to eQTL Catalogue.  
#' If you wish to add custom names to the loci, simply add these as the names of the path list
#'  (e.g. \code{c(BST1="<path>/<to>/<BST1_file>", LRRK2="<path>/<to>/<LRRK2_file>")}). 
#'  Otherwise, loci will automatically named based on their min/max genomic coordinates. 
#'  
#' The minimum columns in these files required to make queries include: 
#' \describe{
#' \item{SNP}{RSID of each SNP.} 
#' \item{CHR}{Chromosome (can be in "chr12" or "12" format).}  
#' \item{POS}{Genomic position of each SNP.}  
#' \item{...}{Optional extra columns.}  
#' }  
#' @param output_dir The folder you want the merged gwas/qtl results to be saved to 
#' (set \code{output_dir=F} if you don't want to save the results).  
#' If \code{split_files=F}, all query results will be merged into one and saved as \emph{<output_dir>/eQTL_Catalogue.tsv.gz}.  
#' If \code{split_files=T}, all query results will instead be split into smaller files and stored in \emph{<output_dir>/}.   
#' @param qtl_search This function will automatically search for any datasets that match your criterion.
#' For example, if you search "Alasoo_2018", it will query the datasets:  
#' \itemize{
#' \item{Alasoo_2018.macrophage_naive}  
#' \item{Alasoo_2018.macrophage_Salmonella}  
#' \item{Alasoo_2018.macrophage_IFNg+Salmonella}
#' } 
#' You can be more specific about which datasets you want to include, 
#' for example by searching: "Alasoo_2018.macrophage_IFNg".
#' You can even search by tissue or condition type (e.g. \code{c("blood","brain")}) 
#' and any QTL datasets containing those substrings (case-insensitive) in their name or metadata will be queried too.  
#' @param use_tabix Tabix is about ~17x faster (\emph{default:} =T) than the REST API (\emph{=F}).  
#' @param nThread The number of CPU cores you want to use to speed up your queries through parallelization.     
#' @param split_files Save the results as one file per QTL dataset (with all loci within each file).  
#' If this is set to \code{=T}, then this function will return the list of paths where these files were saved.
#' A helper function is provided to import and merge them back together in R.  
#' If this is set to  \code{=F}, then this function will instead return one big merged data.table  
#' containing results from all QTL datasets and all loci.  
#' \code{=F} is not recommended when you have many large loci and/or many QTL datasets,  
#' because you can only fit so much data into memory.  
#' @param quant_method eQTL Catalogue actually contains more than just eQTL data.  
#' For each dataset, the following kinds of QTLs can be queried:  
#' \describe{
#' \item{gene expression QTL}{\code{quant_method="ge"} (\emph{default}) or \code{quant_method="microarray"}, depending on the dataset. \strong{catalogueR} will automatically select whichever option is available.}  
#' \item{exon expression QTL}{\emph{*under construction*}  \code{quant_method="ex"}}
#' \item{transcript usage QTL}{\emph{*under construction*}  \code{quant_method="tx"}}
#' \item{promoter, splice junction and 3' end usage QTL}{\emph{*under construction*}  \code{quant_method="txrev"}}   
#' }
#' @param merge_with_gwas Whether you want to merge your QTL query results with your GWAS data 
#' (convenient, but takes up more storage).
#' @param force_new_subset By default, \strong{catalogueR} will use any pre-existing files that match your query. 
#' Set \code{force_new_subset=T} to override this and force a new query.  
#' @param genome_build The genome build of your query coordinates (e.g. \code{gwas_data}). 
#' If your coordinates are in \emph{hg19}, \strong{catalogueR} will automatically lift them over 
#' to \emph{hg38} (as this is the build that eQTL Catalogue uses).  
#' @param progress_bar \code{progress_bar=T} allows progress to be monitored even when multithreading enabled.
#' Requires R package \code{\link{pbmcapply}}.  
#' @param verbose Show more (\code{=T}) or fewer (\code{=F}) messages.  
#' @family eQTL Catalogue 
#' @examples 
#' sumstats_paths <- example_sumstats_paths()
#' 
#' # Merged results
#' # GWAS.QTL <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=1, force_new_subset=T, merge_with_gwas=F, progress_bar=T, split_files=F)
#' # Merged results (parallel)
#' GWAS.QTL <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=4, force_new_subset=T, merge_with_gwas=F, progress_bar=T, split_files=F)
#' 
#' # Split results
#' # gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=1, force_new_subset=T, merge_with_gwas=F, progress_bar=T) 
#' # Split results (parallel)
#' gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths, qtl_search="Alasoo_2018", nThread=4, force_new_subset=T, merge_with_gwas=F, progress_bar=T)
#' GWAS.QTL <- gather_files(file_paths = gwas.qtl_paths)
#' 
#' # Nalls et al example
#' \dontrun{
#' sumstats_paths_Nalls <- list.files("Fine_Mapping/Data/GWAS/Nalls23andMe_2019","Multi-finemap_results.txt", recursive = T, full.names = T)
#' names(sumstats_paths_Nalls) <- basename(dirname(dirname(sumstats_paths_Nalls)))
#' gwas.qtl_paths <- eQTL_Catalogue.query(sumstats_paths=sumstats_paths_Nalls, output_dir="catalogueR_queries/Nalls23andMe_2019", merge_with_gwas=T, nThread=1, force_new_subset=T)
#' }
eQTL_Catalogue.query <- function(sumstats_paths=NULL,
                                 output_dir="./catalogueR_queries",
                                 qtl_search=NULL,
                                 use_tabix=T,
                                 nThread=4, 
                                 # multithread_qtl=T,
                                 # multithread_loci=F,
                                 # multithread_tabix=F,
                                 quant_method="ge",
                                 infer_region=T, 
                                 split_files=T,
                                 merge_with_gwas=T,
                                 force_new_subset=F, 
                                 genome_build="hg19", 
                                 progress_bar=T,
                                 verbose=T){
  # sumstats_paths <-  example_sumstats_paths()
  # merge_with_gwas=F; nThread=4; loci_names = basename(dirname(dirname(sumstats_paths))); qtl.id=qtl_datasets; multithread_qtl=T;  multithread_loci=F; quant_method="ge";   split_files=T;
  # qtl_search =c("ROSMAP","Alasoo_2018","Fairfax_2014", "Nedelec_2016","BLUEPRINT","HipSci.iPSC", "Lepik_2017","BrainSeq","TwinsUK","Schmiedel_2018", "blood","brain")
  # output_dir = "/Volumes/Scizor/eQTL_Catalogue/Nalls23andMe_2019"; qtl.id="Alasoo_2018.macrophage_naive"; force_new_subset=F; progress_bar=F; verbose=T; genome_build="hg19";

  # Aargs hidden from user in favor of automatic optimizer
  # @param multithread_qtl Multi-thread across QTL datasets (good when you're querying lots of QTL datasets).
  # @param multithread_loci Multi-thread across loci (good when you have lots of gwas loci).
  # @param multithread_tabix Multi-thread across within a single tabix file query (good when you have one-several large loci).
  
  query.start <- Sys.time()
  # Setup 
  cleanup_tbi(DIR=dirname(output_dir))  
  lapply_func <- progress_bar_check(progress_bar=progress_bar, 
                                    verbose=verbose)
  # Search metadata for matching datasets 
  qtl_datasets <- eQTL_Catalogue.search_metadata(qtl_search = qtl_search, 
                                                 verbose = verbose)
  
  # Determine multi-threading level
  # multithread_opts <- multithread_handler(multithread_qtl=multithread_qtl, multithread_loci=multithread_loci, multithread_tabix=multithread_tabix)
  multithread_opts <- multithread_optimizer(qtl_datasets=qtl_datasets, 
                                            sumstats_paths=sumstats_paths)
  multithread_qtl<-multithread_opts$qtl; multithread_loci<-multithread_opts$loci; multithread_tabix<-multithread_opts$tabix; 
  
 
  # QUERY eQTL Catalogue
  printer("eQTL_Catalogue:: Querying",length(qtl_datasets),"QTL datasets x",length(sumstats_paths),"GWAS loci",
          paste0("(",length(qtl_datasets)*length(sumstats_paths)," total)"), v=verbose)   
  # ---- Iterate over QTL datasets  
  GWAS.QTL_all <- lapply_func(qtl_datasets, function(qtl_id,
                                                     .sumstats_paths=sumstats_paths,
                                                     .output_dir=output_dir,
                                                     .quant_method=quant_method,
                                                     .infer_region=infer_region,
                                                     .use_tabix=use_tabix,
                                                     .nThread=nThread,
                                                     .split_files=split_files,
                                                     .merge_with_gwas=merge_with_gwas,
                                                     .force_new_subset=force_new_subset,
                                                     .genome_build=genome_build,
                                                     .multithread_loci=multithread_loci,
                                                     .multithread_tabix=multithread_tabix,
                                                     .verbose=verbose){  
    message(qtl_id)
    GWAS.QTL <- NULL
    # try({
      GWAS.QTL <- eQTL_Catalogue.iterate_fetch(sumstats_paths=.sumstats_paths, 
                                               output_dir=.output_dir,
                                               qtl_id=qtl_id,
                                               quant_method=.quant_method,
                                               infer_region=.infer_region, 
                                               use_tabix=.use_tabix, 
                                               nThread=.nThread, 
                                               split_files=.split_files,
                                               merge_with_gwas=.merge_with_gwas, 
                                               force_new_subset=.force_new_subset, 
                                               genome_build=.genome_build,
                                               multithread_loci=.multithread_loci,
                                               multithread_tabix=.multithread_tabix,
                                               progress_bar = F,
                                               verbose=.verbose) 
    # })
   
    return(GWAS.QTL)
  }, mc.cores = if(multithread_qtl) nThread else 1) # END ITERATION OVER QTL_IDS
  
  # Gather paths/results
  if(split_files){
    # Split protocol
    printer("++ Returning list of split files paths.", v=verbose)
    GWAS.QTL_all <- unlist(GWAS.QTL_all) 
  } else{ 
    printer("++ Post-processing merged results.", v=verbose)
    # Merged protocol
    GWAS.QTL_all <- data.table::rbindlist(GWAS.QTL_all, fill = T)
    # Remove NAs 
    GWAS.QTL_all <- subset(GWAS.QTL_all, !is.na(beta.QTL))
    # Check for completeness
    if(!is.null(names(sumstats_paths))){
      missing_queries <- check_missing_queries(qtl_datasets = qtl_datasets, 
                                               loci = names(sumstats_paths), 
                                               GWAS.QTL = GWAS.QTL_all, 
                                               verbose = verbose)
    } 
    # Save
    if(output_dir!=F){
      output_file <- file.path(output_dir,"eQTL_Catalogue.tsv.gz")
      printer("+ Saving merged query results ==>",output_file, v=verbose)
      dir.create(output_dir, showWarnings = F, recursive = T)
      data.table::fwrite(GWAS.QTL_all,
                         file = output_file,
                         nThread = nThread) 
    }  
  } # END MERGED PROTOCOL
  
  # Clean
  cleanup_tbi(DIR=dirname(output_dir))
  # Clock it
  query.end <- Sys.time()
  print(round(query.end-query.start,1))
  return(GWAS.QTL_all)
}  



