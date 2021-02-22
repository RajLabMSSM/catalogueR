


#' Paths to example summary stats
#' 
#' Returns the paths to summary stats stored within \emph{catalogueR}.
#' Each file is the output of a locus that has been fine-mapping using \emph{echolocatoR}.
#' Data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#' 
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @param Rlib_path This function will automatically find your Rlib path, 
#' but you can override this by supplying it manually.
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' These example files can be used 
#' @examples 
#' sumstats_paths <- example_sumstats_paths()
example_sumstats_paths <- function(Rlib_path=NULL){
  if(is.null(Rlib_path)){
    cat_dir <- system.file("extdata/Nalls23andMe_2019",package = "catalogueR")
  } else {
    cat_dir <- file.path(Rlib_path,"catalogueR/extdata/Nalls23andMe_2019")
  } 
  sumstats_paths <- list.files(cat_dir, pattern = "*_subset.tsv.gz", recursive = T, full.names = T)
  locus_names <- unlist(lapply(strsplit(basename(sumstats_paths),"_"),function(x){x[1]}))
  names(sumstats_paths) <- locus_names
  return(sumstats_paths)
}





#' Paths to example eQTL Catalogue query results
#' 
#' Returns the paths to eQTL Catalogue query results stored within \emph{catalogueR}.
#' Each file is a merged data.table of the GWAS summary stats used to make the query, 
#' and the eQTL Catalogue query results (which can contain data for multiple eGenes).
#' 
#' GWAS data originally comes from the Parkinson's disease GWAS
#' by \href{https://www.biorxiv.org/content/10.1101/388165v3}{Nalls et al. (bioRxiv)}.
#' 
#' \describe{
#'   \item{SNP}{SNP RSID}
#'   \item{CHR}{Chromosome}
#'   \item{POS}{Genomic positiion (in basepairs)}
#'   ...
#' }
#' @param Rlib_path This function will automatically find your Rlib path, 
#' but you can override this by supplying it manually.
#' @source \url{https://www.biorxiv.org/content/10.1101/388165v3}
#' @family Nalls23andMe_2019
#' These example files can be used 
#' @examples 
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
example_eQTL_Catalogue_query_paths <- function(Rlib_path=NULL){
  if(is.null(Rlib_path)){
    cat_dir <- system.file("extdata/eQTL_Catalogue_queries",package = "catalogueR")
  } else {
    cat_dir <- file.path(Rlib_path,"catalogueR/extdata/eQTL_Catalogue_queries")
  } 
  sumstats_paths <- list.files(cat_dir, pattern = "*.tsv.gz", recursive = T, full.names = T)
  locus_names <- unlist(lapply(strsplit(basename(sumstats_paths),"_"),function(x){x[1]}))
  names(sumstats_paths) <- locus_names
  return(sumstats_paths)
}






# Consistent means of making split path
make_split_path <- function(output_dir,
                            qtl_id,
                            loc){
  split_path <- file.path(output_dir, 
                          qtl_id, 
                          paste0(loc,"__",qtl_id,".tsv.gz"))
  names(split_path) <- paste(loc,qtl_id,sep="__")
  dir.create(dirname(split_path), showWarnings = F, recursive = T) 
  return(split_path)
}





parse_gwas.qtl_path <- function(gwas.qtl_path,
                                get_locus=F,
                                get_qtl_id=T,
                                sep='__'){
  split_name <- strsplit(basename(gwas.qtl_path),sep)[[1]]
  if(get_locus & get_qtl_id==F)return(split_name[1])
  if(get_locus==F & get_qtl_id)return(gsub(".tsv|.gz","",split_name[2]))
  if(get_locus & get_qtl_id)return(list(locus=split_name[1],
                                        qtl_id=gsub(".tsv|.gz","",split_name[2])))
}




# Interactive data table
createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}




# Interactive data table (when manually creating rmarkdown chunks)
createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}




# rbind a list of files
rbind.file.list <- function(file.list, 
                            verbose=T, 
                            nCores=4){
  merged.dat <- parallel::mclapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x, nThread = 1)
    return(dat)
  }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}




# Print
printer <- function(..., v=T){if(v){print(paste(...))}}





# Make up a locus name based on coordinates
construct_locus_name <- function(gwas_data, 
                                 verbose=T){ 
  printer("++ Constructing locus name from coordinates", v=verbose)
  paste0("locus_chr",gwas_data$CHR[1],"-",min(gwas_data$POS),"-",max(gwas_data$POS)) 
} 





# Try selecting the quant_method requested, but if not available select another
choose_quant_method <- function(ui, 
                                qm, 
                                verbose=T){
  meta <- eQTL_Catalogue.list_datasets(verbose = F)
  meta.sub <- data.frame(subset(meta, unique_id==ui))
  if(qm %in% unique(meta.sub$quant_method) ){
    meta.sub <- subset(meta.sub, quant_method==qm)
  } else {
    meta.sub <- meta.sub[1,]
    printer("+ Selecting quant_method:",meta.sub$quant_method[1], v=verbose)
  }
  return(meta.sub)
} 




# Remove extra .tbi files
cleanup_tbi <- function(DIR="./", 
                        recursive=F){
  tbi_files <- list.files(DIR, pattern = ".gz.tbi$", recursive = recursive)
  if(length(tbi_files)>0)out <- file.remove(tbi_files)
} 




# Get the size of the largest file in a dir
max_file_size <- function(top_dir, recursive=T){
  # 
  files <- list.files(top_dir, full.names = T, recursive = recursive)
  inf <- file.info(files) 
  print(paste("File size info (Mb) for",length(files),"files."))
  summary(inf$size/1e6)
}




# Time a process
timeit <- function(func, digits=3){
  start = Sys.time()
  out <- func
  end = Sys.time()
  # print(paste("Completed in",round(end-start,digits = digits),"seconds."))
  round(end-start,digits = digits)
  return(out)
}




# Add eGene column
add_eGene_col <- function(qtl.dat){
  # Remove any old col
  if("eGene" %in% colnames(qtl.dat)){qtl.dat <- subset(qtl.dat, select=-eGene)}
  qtl.dat <- qtl.dat %>% dplyr::mutate(eGene = ifelse(!is.na(gene.QTL) & gene.QTL!="",gene.QTL, molecular_trait_id.QTL))
  return(qtl.dat)
}




#' Lift genome across builds
#' 
#' @param build_conversion "hg19.to.hg38" (\emph{default}) or "hg38.to.hg19.
#' @family utils
#' @examples 
#' data("BST1")
#' gr.lifted <- liftover(gwas_data=BST1, build.conversion="hg19.to.hg38")
liftover <- function(gwas_data, 
                     build.conversion="hg19.to.hg38",
                     verbose=T){  
  printer("XGR:: Lifting genome build:", build.conversion, v = verbose)
  # Save original coordinates and SNP IDs
  gwas_data <- gwas_data %>% dplyr::mutate(chrom=paste0("chr",gsub("chr","",CHR)),
                                           POS.orig=POS,
                                           SNP.orig=SNP)
  # chain <- rtracklayer::import.chain(con = chain_paths$hg19_to_hg38)
  gr.gwas <- GenomicRanges::makeGRangesFromDataFrame(df =gwas_data, 
                                                     keep.extra.columns = T, 
                                                     seqnames.field = "chrom", 
                                                     start.field = "POS", 
                                                     end.field = "POS") 
  gr.lifted <- XGR::xLiftOver(data.file = gr.gwas, #dplyr::select(gwas_data, CHR, POS), 
                              format.file = "GRanges",
                              build.conversion = build.conversion, 
                              verbose = verbose , 
                              merged = F)  # merge must =F in order to work
  return(gr.lifted)
}




#' Convert HGNC gene symbols to ENSEMBL IDs
#' 
#' @family utils
#' @examples 
#' gene_symbols <- c("BDNF","FOXP2","BST1")
#' ensembl_ids <- hgnc_to_ensembl(gene_symbols)
hgnc_to_ensembl <- function(gene_symbols,
                            unique_only=T,
                            verbose=T){
  printer("++ Converting: HGNC gene symbols ==> Ensembl IDs", v=verbose)
  if(unique_only) gene_symbols <- unique(gene_symbols)
  # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  gene_symbols[is.na(gene_symbols)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = gene_symbols,
                                      keytype = "SYMBOL",
                                      column = "GENEID")
  return(conversion)
}




#' Convert ENSEMBL IDs to HGNC gene symbols
#' 
#' @family utils
#' @examples 
#' ensembl_ids <- c("ENSG00000176697","ENSG00000128573","ENSG00000109743")
#' gene_symbols <- ensembl_to_hgnc(ensembl_ids=ensembl_ids)
ensembl_to_hgnc <- function(ensembl_ids, 
                            unique_only=T,
                            verbose=T){
  printer("++ Converting: Ensembl IDs ==> HGNC gene symbols", v=verbose)
  if(unique_only) ensembl_ids <- unique(ensembl_ids)
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = ensembl_ids,
                                      keytype = "GENEID",
                                      column = "SYMBOL")
  return(conversion)
}




#' Find Consensus SNPs in \emph{echolocatoR} output
#' 
#' @family echolocatoR
#' @examples 
#' data("BST1")
#' BST1 <- find_consensus_SNPs(finemap_dat=BST1)
find_consensus_SNPs <- function(finemap_dat,
                                verbose=T,
                                credset_thresh=.95,
                                consensus_thresh=2,
                                sort_by_support=T,
                                exclude_methods=NULL){
  printer("+ Identifying Consensus SNPs...",v=verbose)
  exclude_methods <- append(exclude_methods,"mean")
  # Find SNPs that are in the credible set for all fine-mapping tools
  CS_cols <- colnames(finemap_dat)[endsWith(colnames(finemap_dat),".CS")]
  CS_cols <- CS_cols[!(CS_cols %in% paste0(exclude_methods,".CS"))]
  if(consensus_thresh=="all"){consensus_thresh<-length(CS_cols)}
  printer("++ support_thresh =",consensus_thresh)
  # Get the number of tools supporting each SNP
  ## Make sure each CS is set to 1
  support_sub <- subset(finemap_dat, select = CS_cols) %>% data.frame()
  support_sub[sapply(support_sub, function(e){e>1})] <- 1
  finemap_dat$Support <- rowSums(support_sub, na.rm = T)
  finemap_dat$Consensus_SNP <- finemap_dat$Support >= consensus_thresh
  # Sort
  if(sort_by_support){
    finemap_dat <- finemap_dat %>% arrange(desc(Consensus_SNP), desc(Support))
  }
  
  # Calculate mean PP
  printer("+ Calculating mean Posterior Probability (mean.PP)...")
  PP.cols <- grep(".PP",colnames(finemap_dat), value = T)
  PP.cols <- PP.cols[!(PP.cols %in% paste0(exclude_methods,".PP"))]
  PP.sub <- subset(finemap_dat, select=c("SNP",PP.cols)) %>% data.frame()# %>% unique()
  PP.sub[is.na(PP.sub)] <- 0
  if(NCOL(PP.sub[,-1]) > 1){
    finemap_dat$mean.PP <- rowMeans(PP.sub[,-1])
  } else{
    finemap_dat$mean.PP <- PP.sub[,-1]
  }
  finemap_dat$mean.CS <- ifelse(finemap_dat$mean.PP>=credset_thresh,1,0)
  
  # PP.sub %>% arrange(desc(mean.PP)) %>% head()
  printer("++",length(CS_cols),"fine-mapping methods used.")
  printer("++",dim(subset(finemap_dat,Support>0))[1],"Credible Set SNPs identified.")
  printer("++",dim(subset(finemap_dat,Consensus_SNP==T))[1],"Consensus SNPs identified.")
  return(finemap_dat)
}



# Make sure there's some way to specify query coordinates
check_coord_input <- function(gwas_data, 
                              chrom,
                              bp_lower,
                              bp_upper){
  if(is.null(gwas_data) & all(c(is.null(chrom), is.null(bp_lower), is.null(bp_upper)))){
    stop("+ User must specify coordinates to fetch using either `gwas_data` OR `chrom`,`bp_lower`, and `bp_upper`")
  }
}

 


# Make sure there's not conflicting levels of parallelization
multithread_handler <- function(multithread_qtl,
                                multithread_loci,
                                multithread_tabix){
  if(sum(c(multithread_qtl,multithread_loci,multithread_tabix))>1){
    warning("++ Only one multithreading option can be used at once. \n", 
            "Setting: `multithread_qtl=T`, `multithread_loci=F`, `multithread_tabix=F`") 
    multithread_qtl=T; multithread_loci <- F; multithread_tabix=F;
  } 
  multithread_opts <- list(qtl=multithread_qtl, 
                           loci=multithread_loci,
                           tabix=multithread_tabix)
  return(multithread_opts)
}




multithread_optimizer <- function(qtl_datasets,
                                  sumstats_paths, 
                                  verbose=T){
  printer("+ Optimizing multi-threading...", v=verbose)
  if(length(qtl_datasets) > length(sumstats_paths)){
    printer("++ Multi-threading across QTL datasets.", v=verbose)
    multithread_opts <- list(qtl=T, 
                             loci=F,
                             tabix=F)
    
  }else {
    printer("++ Multi-threading across loci.", v=verbose)
    multithread_opts <- list(qtl=F, 
                             loci=T,
                             tabix=F)
  }
  return(multithread_opts)
}




update_CS_cols <- function(finemap_dat){
  colnames(finemap_dat) <- gsub("*.Credible_Set$",".CS",colnames(finemap_dat))
  return(finemap_dat)
}




progress_bar_check <- function(progress_bar,
                               verbose=T){
  if(progress_bar){
    if("pbmcapply" %in% installed.packages()){
      lapply_func <- pbmcapply::pbmclapply
    } else {
      printer("++ R package `pbmcapply` not installed. Turning off `progress_bar`.",v=verbose)
      lapply_func <- parallel::mclapply
    }
  } else {
    lapply_func <- parallel::mclapply
  } 
  return(lapply_func)
}



#' Merge files from a list of paths
#' 
#' Merge a list of files into one by stacking them on top of each other (i.e. rbind).
#' @family utils
#' @examples 
#' sumstats_paths <- example_sumstats_paths()
#' merged_dat <- gather_files(file_paths=sumstats_paths)
gather_files <- function(file_paths,
                         nThread=4,
                         verbose=T){
  printer("+ Merging",length(file_paths),"files.", v=verbose)
  printer("+ Using",nThread,paste0(ifelse(nThread>1,"cores.","core.")), v=verbose)
  DAT <- parallel::mclapply(file_paths, function(x){
    dat <- data.table::data.table()
    try({
      if(endsWith(tolower(x),".rds")){
        dat <- readRDS(x)
      } else {
        dat <- data.table::fread(x)
      } 
      dat$file <- basename(x)
    })
  return(dat)
  }, mc.cores = nThread) %>% data.table::rbindlist(fill=T) 
  printer("+ Merged data.table:",nrow(DAT),"rows x",ncol(DAT),"columns.", v=verbose)
  return(DAT)
}





check_missing_queries <- function(qtl_datasets,
                                  loci,
                                  GWAS.QTL,
                                  verbose=T){
  printer("+ Checking results for missing queries...")
  GWASlocus__QTLid.queried <- unlist(lapply(qtl_datasets, function(x){paste(loci,x, sep="__")}))
  GWASlocus__QTLid.found <- unique(paste(GWAS.QTL$Locus.GWAS, GWAS.QTL$qtl_id,sep="__"))
  
  missing_queries <- GWASlocus__QTLid.queried[!GWASlocus__QTLid.queried %in% GWASlocus__QTLid.found]
  if(length(missing_queries)>0){
    printer("+ QTL datasets with no hits/failed to pull data:", v=verbose)
    for(x in missing_queries){printer("  +",x)}
  }  
  return(missing_queries)
}





#' Merge coloc results with SNP-wise GWAS/QTL data
#' 
#' @family coloc
#' @examples 
#' \dontrun{
#' data("coloc_QTLs")
#' query_paths <- example_eQTL_Catalogue_query_paths()
#' output_dir <- unique(dirname(dirname(query_paths))) 
#' } 
retrieve_sumstats_info <- function(output_dir="./catalogueR_queries",
                                   coloc_QTLs){ 
  query_paths <- file.path(output_dir,qtl_id, )
  lapply(1:nrow(coloc_QTLs), function(snp,
                                      .output_dir=output_dir){
    ROW <- coloc_QTLs[i,]
    dat_path <- make_split_path(output_dir = .output_dir,
                                qtl_id = ROW$qtl_id, 
                                loc = ROW$Locus.GWAS)
    dat <- data.table::fread(dat_path, nThread = 1)
  })
  file.path()
  
}



reconstruct_file_names <- function(coloc_QTLs){
  gwas.qtl_files <- unique(paste0(coloc_QTLs$Locus.GWAS,"__",coloc_QTLs$qtl_id,".tsv.gz")) 
  return(gwas.qtl_files)
}


check_dim <- function(df){
  try({
    printer("Data dimensions:",nrow(df),"x",ncol(df))
  }) 
}


#' Infer (effective) sample size from summary stats
#'
#' @family general
#' @keywords internal
#' @examples
#' data("BST1")
#' BST1 <- finemap_DT
#' subset_DT <- get_sample_size(subset_DT = finemap_DT)
get_sample_size <- function(subset_DT,
                            sample_size=NULL,
                            effective_ss=T,
                            verbose=T){
  printer("+ Preparing sample_size (N) column",v=verbose)
  if(!"N" %in% colnames(subset_DT)){
    if(is.null(sample_size)){
      if("N_cases" %in% colnames(subset_DT) & "N_controls" %in% colnames(subset_DT)){
        if(effective_ss){
          printer("++ Computing effective sample size.",v=verbose)
          subset_DT$N <- round(4.0 / (1.0/subset_DT$N_cases + 1.0/subset_DT$N_controls), digits = 0)
        }else {
          printer("++ Inferring sample size from max(N_cases) + max(N_controls):",sample_size,v=verbose)
          subset_DT$N <- subset_DT$N_cases + subset_DT$N_controls
        }
      } else {
        subset_DT$N <- NULL
        printer("++ `sample_size` not provided.",v=verbose)
      }
    } else {
      printer(paste0("++ Using `sample_size = ",sample_size,"` ",v=verbose) );
      subset_DT$N <- sample_size
    }
  } else {printer("+ `N` column already present in data.",v=verbose)}
  return(subset_DT)
}




generate_token <- function(n = 1, len=5) {
  a <- do.call(paste0, replicate(len, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
