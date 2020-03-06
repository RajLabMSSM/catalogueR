### ********************************* ###
### ********************************* ###
##### ------ eQTLcatalogueR ------- #####
### ********************************* ###
### ********************************* ###

# The following functions are derived from the following eQTl Catalogue tutorial
# http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html
# GitHub: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources
# Table of all tabix-indexed QTL datasets:  https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv
# In-depth API documentation: https://www.ebi.ac.uk/eqtl/api-docs/
# Tabix Instructions: https://www.ebi.ac.uk/eqtl/Data_access/
# FTP server: ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv
# Table of tabix files on ftp server: 

# Notes on parallelization
## There's multiple levlels to parallelize on.
## 1. QTL dataset: 
### Speedup
## 2. per Locus query
## 3. Reading in files (e.g. tabix subset directly into R from eQTL Catalogue server)

# You can also get a speedup by using tabix instead of the rest API
## Test: For 3 loci, and X QTL datasets:
### Tabix: 27 seconds (*clear winner! ~17x speedup)
### REST API: 7.5 minutes

# First things first: Loading needed libraries
library("dplyr")
library("ggplot2")
library("readr")
library("stringr")
library("httr")
library("jsonlite")
library("tidyverse")
library("coloc")
library("biomaRt")
library("wiggleplotr")
library("GenomicRanges")
library("biomaRt")
library(AnnotationDbi)
library(data.table)
library(EnsDb.Hsapiens.v75)
library(XGR)

  

# 1. List available eQTL data
catalogueR.list_eQTL_datasets <- function(save_path="./resources",
                                         force_new=F,
                                         verbose=F){
  meta.path <- file.path(save_path,"eQTLcatalogue_tabix_ftp_paths.tsv")
  if(file.exists(meta.path) & force_new==F){
    printer("+ Importing saved metadata.",v = verbose)
    meta <- data.table::fread(meta.path, nThread = 4) 
    subset(meta, unique_id=="Schmiedel_2018.Treg_naive") 
    meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path)) 
  } else {
    printer("+ Downloading metadata from server.",v = verbose)
    URL <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
    meta <- data.table::fread(URL, nThread = 4)
    meta <- meta %>% dplyr::transmute(unique_id=paste0(study,".",qtl_group), !!!.) 
    meta <- meta %>% dplyr::mutate(ftp_path= gsub("Fairfax_2014_monocyte","Fairfax_2014",ftp_path)) 
    if(save_path!=F){
      printer("Saving metadata ==>",meta.path)
      if(!dir.exists(dirname(meta.path)))dir.create(dirname(meta.path))
      data.table::fwrite(meta, meta.path, sep="\t")
    }
  }
  printer(length(unique(meta$study)),"unique eQTL datasets:",v=verbose)
  # printer(unique(meta$study))
  return(meta)
}


 
# 2. Get QTL data by region

# 2.1 Method 1: Tabix
## Faster alternative to REST API
catalogueR.fetch_tabix <- function(unique_id,
                                      quant_method="ge",
                                      infer_region=T,
                                      gwas_data=NULL,
                                      chrom=NULL,
                                      bp_lower=NULL,
                                      bp_upper=NULL,
                                      is_gwas=F,
                                      nThread=4){
  # quant_method="ge"; infer_region=T;is_gwas=F; remove_tmp=F;  add_chr=T
  tabix.start = Sys.time()
  # Get region
  if(infer_region & !is.null(gwas_data)){
    print("+ Inferring coordinates from gwas_data")
    chrom <- unique(gwas_data$CHR)
    if(length(chrom)>1){stop("More than one chromosome detected.")}
    bp_lower <- min(gwas_data$POS)
    bp_upper <- max(gwas_data$POS)
  }
  region <- paste0(chrom,":",bp_lower,"-",bp_upper)
  # printer("+ TABIX:: Querying region:", region)
  # Get metadata
  # Rename var to avoid issues with subsetting
  ui <-unique_id
  qm <- quant_method
  meta <- catalogueR.list_eQTL_datasets(force_new = F, save_path = F)
  meta.sub <- subset(meta, unique_id==ui) %>% data.frame()
  if(qm %in% unique(meta.sub$quant_method) ){
    meta.sub <- subset(meta.sub, quant_method==qm)
  } else {
    meta.sub <- subset(meta.sub, quant_method==meta.sub$quant_method[1])
    printer("+ Selecting quant_method:",meta.sub$quant_method[1])
  }
  # Run tabix
  ## Get the header (according to email with Kaur Alassoo)
  catalogueR.tabix_header <- function(tabix_path){
    eqtl_cat_path <- "./resources"
    tabix_headers <- file.path(eqtl_cat_path,"tabix_header.txt")
    if(file.exists(tabix_headers)){
      header <- as.list(data.table::fread(tabix_headers))[[1]]
    } else {
      dir.create(eqtl_cat_path,showWarnings = F, recursive = T)
      header.path <- paste("curl -s",tabix_path,"| zcat | head -n 1")
      header <-  colnames(data.table::fread(cmd = header.path))
      data.table::fwrite(list(header), file=tabix_headers)
    }
    return(header)
  }
  header <- catalogueR.tabix_header(tabix_path = meta.sub$ftp_path)

  # Run tabix
  # Read directly into R rather than saving tabix subset
  qtl.subset <- data.table::fread(cmd=paste("tabix",
                                     # "--print-header",
                                     meta.sub$ftp_path,
                                     region),
                           nThread = nThread)
  colnames(qtl.subset) <- c("Locus.QTL",paste0(header,".QTL"))
  tabix.end = Sys.time() 
  printer("eQTL_catalogue::",nrow(qtl.subset),"eSNPs returned in", round(as.numeric(tabix.end-tabix.start),1),"seconds.")  
  return(qtl.subset)
}



# 2.2 Method 2: RESTful API
### Slower than tabix?
fetch_from_eqtl_cat_API <- function(link,
                                    is_gwas = FALSE){
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
  is_paginated <- !str_detect(link,"paginate=False")
  # message("isPagined:", is_paginated)
  page = 1
  merged_summaries <- data.frame()
  while(!is.null(link)){
    # print(paste0("Fetching page #",page))
    api_raw_data <- jsonlite::fromJSON(link, simplifyDataFrame = TRUE, flatten = TRUE)
    link <- api_raw_data$`_links`$`next`$href
    if (is_empty(api_raw_data$`_embedded`$associations)) {
      return(merged_summaries)
    }
    eqtl_raw_list_data <- do.call(rbind, lapply(api_raw_data$`_embedded`$associations, rbind))
    eqtl_data <- nullToNA(eqtl_raw_list_data) %>% as.matrix() %>% as_tibble()
    if (is_paginated) { eqtl_data <- dplyr::select(eqtl_data, -c("_links")) }
    eqtl_data <- tidyr::unnest(eqtl_data, cols = cols_to_nest)
    if (!is.null(link)) {
      page <- page + 1
    }
    merged_summaries <- merged_summaries %>% rbind(eqtl_data)
  }
  return(merged_summaries[cols_to_nest])
}




catalogueR.fetch_restAPI <- function(unique_id, #Alasoo_2018.macrophage_naive
                                     quant_method="ge",
                                     infer_region=T,
                                     gwas_data=NULL,
                                     chrom=NULL,
                                     bp_lower=NULL,
                                     bp_upper=NULL,
                                     is_gwas=F, # refers to the datasets being queried
                                     size=NULL) {
  restAPI.start = Sys.time()
  # Get region
  if(infer_region & !is.null(gwas_data)){
    print("+ Inferring coordinates from gwas_data")
    chrom <- gsub("chr","",unique(gwas_data$CHR)[1])
    bp_lower <- min(gwas_data$POS)
    bp_upper <- max(gwas_data$POS)
  }
  # Get metadata
  ui <- unique_id
  qm <- quant_method
  meta <- catalogueR.list_eQTL_datasets(force_new = F, save_path = F)
  meta.sub <- subset(meta, unique_id==ui) %>% data.frame()
  if(qm %in% unique(meta.sub$quant_method) ){
    meta.sub <- subset(meta.sub, quant_method==qm)
  } else {
    meta.sub <- subset(meta.sub, quant_method==meta.sub$quant_method[1])
    printer("+ Selecting quant_method:",meta.sub$quant_method[1])
  }
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
                 ifelse(!is.null(size),paste0("&size=",size),""))

  message("Fetching: ", link)
  qtl.subset <- fetch_from_eqtl_cat_API(link = link,  is_gwas = is_gwas)
  colnames(qtl.subset) <-  paste0(colnames(qtl.subset),".QTL")  
  restAPI.end = Sys.time()
  printer("eQTL_catalogue::",nrow(qtl.subset),"eSNPs returned in", round(as.numeric(restAPI.end-restAPI.start),1),"seconds.")  
  return(qtl.subset)
}

catalogueR.fetch <- function(unique_id,
                             quant_method="ge",
                             infer_region=T,
                             gwas_data=NULL,
                             is_gwas=F,
                             nThread=1,
                             use_tabix=T,
                             chrom=NULL,
                             bp_upper=NULL,
                             bp_lower=NULL){
  if(use_tabix){
    # Tabix is about ~17x faster than the REST API.
    gwas.qtl <- catalogueR.fetch_tabix(unique_id=unique_id,
                                       quant_method="ge",
                                       infer_region=T,
                                       gwas_data=gwas_data,
                                       is_gwas=F,
                                       nThread = 1,)
  } else {
    gwas.qtl <- catalogueR.fetch_restAPI(unique_id=unique_id,
                                         quant_method="ge",
                                         infer_region=T,
                                         gwas_data=gwas_data,
                                         is_gwas=F)
  }
  return(gwas.qtl)
}

# 4. Run coloc
# Method to perform colocalisation analysis. 



# 5. Merge results
make_assocs_list_to_merged_plottable <- function(all_coloc_dt, quant_method = "ge"){
  merged_assoc_data <- data.frame()
  for (index in 1:(all_coloc_dt$eqtl_assocs_in_region %>% length())) {
    assoc_df <- all_coloc_dt$eqtl_assocs_in_region[[index]]
    assoc_df$coloc_PP4 <- round(all_coloc_dt$colocs[6,index], 3)
    assoc_df$coloc_PP3 <- round(all_coloc_dt$colocs[5,index], 3)
    assoc_df$track <- paste0("eQTL Catalogue\n",
                             assoc_df$study_id, "_",quant_method,"\n",
                             assoc_df$qtl_group,"\n",assoc_df$molecular_trait_id)
    if (index==1) {
      merged_assoc_data <- assoc_df
    } else {
      merged_assoc_data <- merged_assoc_data %>% rbind(assoc_df)
    }
  }
  return(merged_assoc_data)
}


catalogueR.merge_gwas_qtl <- function(gwas_data, 
                                      qtl.subset){ 
  # Merging and allele flipping
  gwas.qtl <- data.table:::merge.data.table(x = data.table::data.table(gwas_data),
                                            y = data.table::data.table(qtl.subset),
                                            # all.x = T,
                                            by.x = c("SNP"), # effect_allele
                                            by.y = c("rsid.QTL") ) %>%
    dplyr::mutate(effect.is.ref=ifelse(A1==ref.QTL,T,F),
                  effect.is.alt=ifelse(A2==alt.QTL,T,F) ) %>%
    # subset(effect.is.ref|effect.is.alt) %>%
    data.table::data.table()
  return(gwas.qtl)
}


catalogueR.eQTL_catalogue.iterate_loci <- function(sumstats_paths, 
                                                   loci_names=NULL,
                                                   output_path,
                                                    qtl_id,
                                                    quant_method="ge",
                                                    infer_region=T, 
                                                    use_tabix=T,
                                                    multithread_loci=T,
                                                    nThread=4,
                                                    split_files=T,
                                                    merge_with_gwas=F,
                                                    force_new_subset=F,
                                                    progress_bar=T,
                                                    genome_build="hg19"){
  apply_func <- ifelse(progress_bar, pbmcapply::pbmclapply, parallel::mclapply)
  GWAS.QTL <-  apply_func(sumstats_paths, function(loc_path){  
    # Import GWAS data
    gwas_data <- data.table::fread(loc_path)
    gwas_data$CHR <- gsub("chr","",gwas_data$CHR)   # get rid of "chr" just in case
    
    # Name Locus
    if("Locus" %in% colnames(gwas_data)){
      printer("++ Using GWAS locus as locus name.")
      loc <- unique(gwas_data$Locus)[1]
    } else {
      if(length(loci_names)==length(sumstats_paths)){
        i <- setNames(1:length(sumstats_paths),sumstats_paths)[[loc_path]]
        loc <- loci_names[i]
      }
      if(length(loci_names)!=length(sumstats_paths)|is.null(loci_names)){
        printer("+ `loci_names` and `sumstats_paths` are of different lengths.",
                "Constructing new loci names instead.")
        loc <- construct_locus_name(gwas_data)
      }
      gwas_data <- cbind(Locus=loc, gwas_data)
    }
    
    message("_+_+_+_+_+_+_+_+_--- Locus ",loc)
    # Test if query file already exists
    split_path <- file.path(output_path, qtl_id, paste0(loc,"_locus__&__",qtl_id,".tsv.gz"))
    dir.create(dirname(split_path), showWarnings = F, recursive = T) 
    
    if(file.exists(split_path) & force_new_subset==F){
      printer("++ Using pre-existing file...")
      qtl.subset <- data.table::fread(split_path)
    } else {
      
      # Convert from GRCh37 to GRCh38
      if(genome_build %in% c('hg19','hg18')){
        gr.lifted <- XGR.liftover(gwas_data, build.conversion = paste0(genome_build,".to.hg38"))
        gwas_data <- data.frame(gr.lifted)  %>% 
          dplyr::mutate(CHR=gsub("chr","",seqnames),
                        POS=start) %>% 
          dplyr::select(-c("seqnames","start","end","width","strand")) %>%
          data.table::as.data.table()
      }
      
      gwas.qtl <- data.table::data.table()
      try({  
        qtl.subset <- catalogueR.fetch(unique_id = qtl_id,
                                       quant_method=quant_method,
                                       infer_region=infer_region,
                                       gwas_data=gwas_data,
                                       is_gwas=F,
                                       nThread=1,
                                       use_tabix=use_tabix,
                                       chrom=NULL,
                                       bp_upper=NULL,
                                       bp_lower=NULL) 
        # Merge results
        if(merge_with_gwas){
          gwas.qtl <- catalogueR.merge_gwas_qtl(gwas_data, qtl.subset)
        } else {gwas.qtl <- qtl.subset}  
        # Add locus name
        gwas.qtl <- cbind(Locus.GWAS=loc, gwas.qtl)
        # Get QTL gene names 
        gene_dict <- ensembl_to_hgnc(ensembl_ids = gwas.qtl$gene_id.QTL)
        gwas.qtl$gene.QTL <- gene_dict[gwas.qtl$molecular_trait_object_id.QTL] 
      })  
      # Save
      if(split_files){data.table::fwrite(gwas.qtl, split_path, sep="\t")}
    }  
    # Return
    if(split_files){return(split_path)} else {return(gwas.qtl)} 
  }, mc.cores = ifelse(multithread_loci,nThread,1))
  # Return
  if(split_files){return(unlist(GWAS.QTL))} else {
    GWAS.QTL <- data.table::data.table(qtl.id=qtl_id, GWAS.QTL) %>% data.table::rbindlist(fill = T)
    return(GWAS.QTL)
  } 
  message(" ")
}



XGR.liftover <- function(gwas_data, 
                         build.conversion="hg19.to.hg38",
                         verbose=F){  
  printer("XGR:: Lifting genome build:", build.conversion, v = verbose)
  # Save original coordinates and SNP IDs
  gwas_data <- gwas_data %>% dplyr::mutate(chrom=paste0("chr",CHR),
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


catalogueR.run <- function(sumstats_paths=NULL,
                           loci_names=NULL,
                           output_path="./example_data/Nalls23andMe_2019",
                           qtl_search=NULL,
                           use_tabix=T,
                           nThread=4, 
                           multithread_qtl=T,
                           multithread_loci=F,
                           quant_method="ge",
                           infer_region=T, 
                           split_files=T,
                           merge_with_gwas=T,
                           force_new_subset=F,
                           progress_bar=T,
                           genome_build="hg19"){
  library("dplyr")
  library("ggplot2")
  library("readr")
  library("stringr")
  library("httr")
  library("jsonlite")
  library("tidyverse")
  library("coloc")
  library("biomaRt")
  library("wiggleplotr")
  library("GenomicRanges")
  library("biomaRt")
  # sumstats_paths <- list.files("./Data/GWAS/Nalls23andMe_2019",  pattern = "*Multi-finemap_results.txt|*Multi-finemap.tsv.gz",  recursive = T, full.names = T) 
  # merge_with_gwas=F; nThread=4; loci_names = basename(dirname(dirname(sumstats_paths))); qtl.id=qtl_datasets; multithread_qtl=T;  multithread_loci=F; quant_method="ge";   split_files=T;
  # qtl_search =c("ROSMAP","Alasoo_2018","Fairfax_2014", "Nedelec_2016","BLUEPRINT","HipSci.iPSC", "Lepik_2017","BrainSeq","TwinsUK","Schmiedel_2018", "blood","brain")
  # output_path = "/Volumes/Scizor/eQTL_catalogue/Nalls23andMe_2019"; qtl.id="Alasoo_2018.macrophage_naive"; force_new_subset=F;
  
  # Helper functions
  construct_locus_name <- function(gwas_data){ 
    paste0("locus_",gwas_data$CHR[1],":",min(gwas_data$POS),"-",max(gwas_data$POS)) 
  } 
  cleanup_tbi <- function(DIR="./"){
    tbi_files <- list.files(DIR, pattern = ".gz.tbi$")
    if(length(tbi_files)>0)file.remove(tbi_files)
  }
  
  # Setup
  cleanup_tbi(DIR="./")
  if(multithread_qtl & multithread_loci){
    printer("++ `multithread_qtl` and `multithread_loci` can't both =TRUE. Setting `multithread_loci=F`");
    multithread_loci <- F
  } 
  # Search metadata for matching datasets
  meta <- catalogueR.list_eQTL_datasets(force_new = F)
  if(is.null(qtl_search)){
    printer("eQTL_catalogue:: Gathering data for all QTL Catalogue datasets...")
    qtl_datasets <- unique(meta$unique_id)
  }else {
    qtl_datasets <- grep(pattern = paste(qtl_search,collapse="|"),
                         x = meta$unique_id,
                         value = T,
                         ignore.case = T) %>% unique()
  }
  
  # QUERY eQTL Catalogue
  printer("eQTL_catalogue:: Querying",length(qtl_datasets),"QTL datasets x",length(sumstats_paths),"GWAS loci.")  
  { start_time <- Sys.time()
    # ---- Iterate over QTL datasets  
    GWAS.QTL_all <- parallel::mclapply(qtl_datasets, function(qtl.id){ 
        # qtl.start = Sys.time()
        message(qtl.id)
        try({ 
          # ---- Iterate over GWAS loci 
          GWAS.QTL <- catalogueR.eQTL_catalogue.iterate_loci(sumstats_paths=sumstats_paths,
                                                             loci_names=loci_names, 
                                                             output_path=output_path,
                                                             qtl_id=qtl.id,
                                                             quant_method=quant_method,
                                                             infer_region=infer_region, 
                                                             use_tabix=use_tabix,
                                                             multithread_loci=multithread_loci,
                                                             nThread=nThread, 
                                                             split_files=split_files,
                                                             merge_with_gwas=merge_with_gwas, 
                                                             force_new_subset=force_new_subset,
                                                             progress_bar=progress_bar,
                                                             genome_build=genome_build)
        })
        # qtl.end = Sys.time()
        # printer("+ Completed queries in",as.numeric(round(qtl.end-qtl.start,1)),"seconds.")
      return(GWAS.QTL)
    }, mc.cores = ifelse(multithread_qtl,nThread,1)) 
    
    # Gather results
    if(split_files){
      printer("++ Returning list of split files paths.")
      return(unlist(GWAS.QTL_all))
    } else{ return(data.table::rbindlist(GWAS.QTL_all, fill = T))  }
    end_time <- Sys.time()
    print(end_time - start_time)
  }

  cleanup_tbi(DIR="./")
  
  # Get QTL gene names 
  # gene_dict <- ensembl_to_hgnc(ensembl_ids = GWAS.QTL_all$gene_id.QTL)
  # GWAS.QTL_all$gene.QTL <- gene_dict[GWAS.QTL_all$molecular_trait_object_id.QTL]
  
  GWAS.QTL_filt <- subset(GWAS.QTL_all, !is.na(beta.QTL))
  missing.qtls <- qtl_datasets[!qtl_datasets %in% unique(GWAS.QTL_filt$qtl.id)]
  if(length(missing.qtls)>0){
    printer("eQTL_Catalogue: QTL datasets with no hits/failed to pull data:")
    for(x in missing.qtls){printer("  +",x)}
  }
  printer("Saving merged query results ==>",output_path)
  if(!dir.exists(dirname(output_path)))dir.create(dirname(output_path))
  data.table::fwrite(GWAS.QTL_filt,
                     file = file.path(output_path,"eQTL_Catalogue.tsv.gz"),
                     nThread = nThread) 
  return(GWAS.QTL_all)
} ### End main function


catalogueR.gather_results <- function(qtl.paths, 
                                      qtl_pval_filter=.05){
  qtl.paths <- list.files("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_catalogue/", 
                          pattern="*.query.tsv.gz$", full.names = T)
  qtl = rbind.file.list(qtl.paths[1:2]) 
}
  


## ---------------- General Functions ----------------  ##

rbind.file.list <- function(file.list, 
                            verbose=T, 
                            nCores=4){
  merged.dat <- parallel::mclapply(file.list, function(x){
    printer(x, v = verbose)
    dat <- data.table::fread(x)
    return(dat)
  }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)
  return(merged.dat)
}

printer <- function(..., v=T){if(v){print(paste(...))}}

hgnc_to_ensembl <- function(gene_symbols){
  # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  gene_symbols[is.na(gene_symbols)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = gene_symbols,
                                      keytype = "SYMBOL",
                                      column = "GENEID")
  return(conversion)
}
ensembl_to_hgnc <- function(ensembl_ids){
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = ensembl_ids,
                                      keytype = "GENEID",
                                      column = "SYMBOL")
  return(conversion)
}

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

createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}




catalogueR.get_colocs <- function(qtl.egene, 
                                  gwas.region,
                                  merge_by_rsid=T){
  # http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html
  # Subset to overlapping SNPs only
  if(merge_by_rsid){
    shared = intersect(qtl.egene$SNP, gwas.region$SNP)
    eqtl_shared = dplyr::filter(qtl.egene, SNP %in% shared) %>% 
      dplyr::mutate(variant_id = SNP) %>% unique()
    gwas_shared = dplyr::filter(gwas.region, SNP %in% shared) %>% 
      dplyr::mutate(variant_id = SNP) %>% unique()
  } else {
    shared = intersect(qtl.egene$position.QTL, gwas.region$POS)
    eqtl_shared = dplyr::filter(qtl.egene, position.QTL %in% shared) %>% 
      dplyr::mutate(variant_id = as.character(position.QTL))
    gwas_shared = dplyr::filter(gwas.region, POS %in% shared) %>% 
      dplyr::mutate(variant_id = as.character(POS))
  }
  
  if(length(shared)==0){
    message("catalogueR:COLOC:: No SNPs shared between GWAS and QTL subsets.")
    coloc_res <- list(summary="No SNPs shared between GWAS and QTL subsets.",
                      results=data.table::data.table(    snp=NA,
                                                         pvalues.df1=NA,
                                                         MAF.df1=NA,
                                                         V.df1=NA,
                                                         z.df1=NA,
                                                         r.df1=NA,
                                                         lABF.df1=NA,
                                                         V.df2=NA,
                                                         z.df2=NA,
                                                         r.df2=NA,
                                                         lABF.df2=NA,
                                                         internal.sum.lABF=NA),
                      Locus=gwas_shared$Locus[1])

  } else { 
    # RUN COLOC
    # QTL data
    eQTL_dataset = list(pvalues = eqtl_shared$pvalue.QTL, 
                        N = (eqtl_shared$an.QTL)[1]/2, #The sample size of the eQTL dataset was 84
                        MAF = eqtl_shared$maf.QTL, 
                        type = "quant", 
                        beta = eqtl_shared$beta.QTL,
                        snp = eqtl_shared$variant_id)
    # GWAS data
    gwas_dataset = list(beta = gwas_shared$Effect, #If log_OR column is full of NAs then use beta column instead
                        varbeta = gwas_shared$StdErr^2, 
                        type = "cc", 
                        snp = gwas_shared$variant_id,
                        s = 0.5, #This is actually not used, because we already specified varbeta above.
                        MAF = gwas_shared$MAF)
    
    coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, 
                                 dataset2 = gwas_dataset,
                                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    coloc_res$Locus <- gwas_shared$Locus[1]
    report <- COLOC.report_summary(coloc.res = coloc_res, 
                                   PP_threshold = .8)
    }
  return(coloc_res)
}




catalogueR.run_coloc <- function(gwas.qtl_paths,
                                 # gwas_paths=NULL,
                                 # qtl_paths=NULL,
                                 nCores=4){
  # qtl.ID=unique(qtl.dat$qtl.id)[1]; eGene=unique(qtl.dataset$gene.QTL)[1];
  # qtl.ID = "Fairfax_2014.monocyte_naive"; eGene = "TSC22D2"; qtl.path=gwas.qtl_paths[1]
  # qtl.path="/Volumes/Scizor/eQTL_catalogue/Nalls23andMe_2019/Alasoo_2018.macrophage_IFNg/MCCC1_locus__&__Alasoo_2018.macrophage_IFNg.tsv.gz"

  # ---- Iterate over QTL datasets 
  coloc_qtls <- lapply(unique(gwas.qtl_paths), function(qtl.path){
    qtl.ID <- strsplit(gsub(".tsv.gz","", basename(qtl.path)), "__&__")[[1]][2]
    printer("+ QTL Dataset =",qtl.ID)
    qtl.dat <- data.table::fread(qtl.path, nThread = 4)
    if(!"qtl.id" %in% colnames(qtl.dat)){qtl.dat <- cbind(qtl.dat, qtl.id=qtl.ID)}
    qtl.dataset <- subset(qtl.dat, qtl.id==qtl.ID & !is.na(gene.QTL) & gene.QTL!="")
    if("Effect" %in% colnames(qtl.dat)){
      gwas.region <- qtl.dataset %>% 
        dplyr::select(Locus,Locus.GWAS, SNP, CHR,POS, P, Effect, StdErr, Freq, MAF, N_cases, N_controls, proportion_cases, A1, A2)
    }  
       
    # ---- Iterate over QTL eGenes
    coloc_eGenes <- parallel::mclapply(unique(qtl.dataset$gene.QTL), function(eGene){ 
      printer("+++ eGene =",eGene)
      qtl.egene <- subset(qtl.dataset, gene.QTL==eGene) 
      # gwas.region <- subset(gwas_data, SNP %in% qtl.egene$rsid.QTL) 
      coloc_res <- catalogueR.get_colocs(qtl.egene = qtl.egene, 
                                         gwas.region = gwas.region, 
                                         merge_by_rsid = T)
      coloc_summary <- as.list(coloc_res$summary)
      
      coloc_DT <- data.table::data.table(Locus.GWAS=coloc_res$Locus, 
                                         qtl.id=qtl.ID,
                                         eGene=eGene,
                                         coloc_res$results,
                                         PP.H0=coloc_summary$PP.H0.abf,
                                         PP.H1=coloc_summary$PP.H1.abf,
                                         PP.H2=coloc_summary$PP.H2.abf,
                                         PP.H3=coloc_summary$PP.H3.abf,
                                         PP.H4=coloc_summary$PP.H4.abf)
      return(coloc_DT)
    }, mc.cores = nCores) %>% data.table::rbindlist(fill=T)  
    return(coloc_eGenes) 
    
  }) %>% data.table::rbindlist(fill=T)
  return(coloc_qtls)
}









# # Method to get significant associations with lead GWAS variant ID.
# get_significant_assocs_of_var <- function(study_ids,
#                                           variant_id,
#                                           p_upper = 0.0001,
#                                           quant_method = "ge"){
#   # In-depth API tutorial: https://www.ebi.ac.uk/eqtl/api-docs/
#   rnaseq_studies_sign_assocs <- data.frame()
#   for(study_name in study_ids){
#     print(study_name)
#     ge_study_query_gwas_lead <- paste0("https://www.ebi.ac.uk/eqtl/api/associations?variant_id=", variant_id,
#                                        "&p_upper=", p_upper,
#                                        "&study=", study_name,
#                                        "&quant_method=", quant_method)
#     fetched_df <- catalogueR.fetch_data(link = ge_study_query_gwas_lead)
#     rnaseq_studies_sign_assocs <- rnaseq_studies_sign_assocs %>% rbind(fetched_df)
#   }
#   if(nrow(rnaseq_studies_sign_assocs)>0){
#     rnaseq_studies_sign_assocs$quant_method <- quant_method
#   } else {print(paste("+ No associations at p <",p_upper,"found."))}
#   return(rnaseq_studies_sign_assocs)
# }
# 
# # lead_var_gwas = dplyr::arrange(gwas_data, p_value)[1,]
# # rnaseq_study_names <- c("Alasoo_2018", "BLUEPRINT", "BrainSeq", "GENCORD", "GEUVADIS",
# #                         "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016",
# #                         "Schwartzentruber_2018", "TwinsUK", "van_de_Bunt_2015")
# # rnaseq_studies_sign_assocs <- get_significant_assocs_of_var(study_ids = rnaseq_study_names,
# #                                                             variant_id = lead_var_gwas$variant_id)
# 
# 
# # Util method for saving plots in different formats
# # save_ggplots <- function(plot, path = ".", filename = "unnamed_plot", height = 15, width = 15){
# #   ggsave(plot = plot,
# #          filename = paste0(filename, ".eps"),
# #          path = path,
# #          device = "eps",
# #          height = height,
# #          width = width,
# #          units = "cm",
# #          dpi = 300)
# #
# #   ggsave(plot = plot,
# #          filename = paste0(filename, ".png"),
# #          path = path,
# #          device = "png",
# #          height = height,
# #          width = width,
# #          units = "cm",
# #          dpi = 300)
# #
# #   ggsave(plot = plot,
# #          filename = paste0(filename, ".pdf"),
# #          path = path,
# #          device = "pdf",
# #          height = height,
# #          width = width,
# #          units = "cm",
# #          dpi = 300)
# # }
# 
# 
# # Method to prepare data for plotting a faceted figure
# make_assocs_list_to_merged_plottable <- function(all_coloc_dt, quant_method = "ge"){
#   merged_assoc_data <- data.frame()
#   for (index in 1:(all_coloc_dt$eqtl_assocs_in_region %>% length())) {
#     assoc_df <- all_coloc_dt$eqtl_assocs_in_region[[index]]
#     assoc_df$coloc_PP4 <- round(all_coloc_dt$colocs[6,index], 3)
#     assoc_df$coloc_PP3 <- round(all_coloc_dt$colocs[5,index], 3)
#     assoc_df$track <- paste0("eQTL Catalogue\n",
#                              assoc_df$study_id, "_",quant_method,"\n",
#                              assoc_df$qtl_group,"\n",assoc_df$molecular_trait_id)
# 
#     if (index==1) {
#       merged_assoc_data <- assoc_df
#     } else {
#       merged_assoc_data <- merged_assoc_data %>% rbind(assoc_df)
#     }
#   }
#   return(merged_assoc_data)
# }
# 
# 
# 
# 
# plot_faceted_multi_manhattans <- function(merged_eqtl_assocs,
#                                           gwas_data_for_trait,
#                                           save_plot = FALSE,
#                                           save_dir = ".",
#                                           save_filename = "new_manhattan_plot",
#                                           save_width = 15,
#                                           save_height = NA,
#                                           no_GWAS_plot = FALSE) {
#   # get shared positions between GWAS and eQTL data
#   shared_positions = intersect(merged_eqtl_assocs$position,
#                                gwas_data_for_trait$base_pair_location)
# 
#   eqtl_shared = dplyr::filter(merged_eqtl_assocs, position %in% shared_positions) %>%
#     dplyr::mutate(variant_id = as.character(position))
#   gwas_shared = dplyr::filter(gwas_data_for_trait, base_pair_location %in% shared_positions) %>%
#     dplyr::mutate(variant_id = as.character(base_pair_location))
# 
#   merged_data_for_plot <- as.data.frame(eqtl_shared %>% dplyr::select(pvalue, position, track))
#   if (!no_GWAS_plot) {
#     gwas_shared$track <- "GWAS RA"
#     gwas_sagred_trans <- gwas_shared %>%
#       dplyr::select(p_value, base_pair_location, track) %>%
#       stats::setNames(c("pvalue", "position", "track"))
#     merged_data_for_plot <- merged_data_for_plot %>% rbind(gwas_sagred_trans)
# 
#     merged_data_for_plot$track <- factor(merged_data_for_plot$track) %>%
#       forcats::fct_relevel("GWAS RA", after = 0)
#   }
# 
#   region_coords = c(45980000, 46200000)
#   plot_base_all = ggplot(merged_data_for_plot, aes_(x = ~ position, y = ~
#                                                       -log(pvalue, 10))) + geom_blank()
# 
#   plot_gwas_eqtl_plot = plot_base_all +
#     geom_point(color = "deepskyblue4") +
#     theme_light() +
#     ylab(expression(paste("-", log[10], " p-value"))) +
#     scale_x_continuous(limits = region_coords, expand = c(0, 0)) +
#     facet_grid(track ~ ., scales = "free_y") +
#     theme(
#       plot.margin = unit(c(0.1, 1, 0.1, 1), "line"),
#       axis.text.x = element_blank(),
#       axis.title.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       legend.position = "none",
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.text.y = element_text(colour = "grey10"),
#       strip.background = element_rect(fill = "grey85")
#     )
# 
#   summary_coloc <- merged_eqtl_assocs %>%
#     group_by(track) %>%
#     summarise(coloc_PP4 = unique(coloc_PP4), coloc_PP3 = unique(coloc_PP3))
# 
#   dat_text <- data.frame(
#     label = c(paste0("PP4: ", summary_coloc$coloc_PP4,
#                      "\nPP3: ", summary_coloc$coloc_PP3)),
#     track = c( summary_coloc$track))
# 
#   if (!no_GWAS_plot) {
#     dat_text <- rbind(data.frame(label="", track="GWAS RA"), dat_text)
#   }
# 
#   plot_gwas_eqtl <- plot_gwas_eqtl_plot + geom_text(
#     data    = dat_text,
#     mapping = aes(x = -Inf, y = -Inf, label = label),
#     hjust   = -0.1,
#     vjust = -3.4
#   )
# 
#   if (save_plot) {
#     save_height= ifelse(test = is.na(save_height),
#                         yes = ((as.integer(!no_GWAS_plot) + unique(merged_eqtl_assocs$track) %>% length()) * 4),
#                         no = save_height)
#     message("Saving: ", save_filename)
#     save_ggplots(
#       plot = plot_gwas_eqtl,
#       path = save_dir,
#       filename = save_filename,
#       height = save_height,
#       width = save_width
#     )
#   }
#   return(plot_gwas_eqtl)
# }
# 
# 
# # Here we start the analysis
# 
# ### We first fetch GWAS data
# # Chromosome: 20
# # study_accession: GCST002318
# # bp_lower: 45980000
# # bp_upper: 46200000
# catalogueR.full_example <- function(){
# 
#   RA_gwas_query_str <- "https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/20/associations?study_accession=GCST002318&bp_lower=45980000&bp_upper=46200000&size=1000"
#   gwas_data <- fetch_from_eqtl_cat_API(link = RA_gwas_query_str, is_gwas = TRUE)
#   message("Downloaded ", gwas_data %>% nrow(), " associations from GWAS Catalogue")
#   gwas_data <- gwas_data %>%
#     dplyr::filter(!ci_lower %>% is.na()) %>%
#     dplyr::filter(!ci_upper %>% is.na()) %>%
#     dplyr::mutate(log_OR = log(odds_ratio)) %>%
#     dplyr::mutate(se = (log(ci_upper)-log(ci_lower))/3.92)
#   message("There remain ", gwas_data %>% nrow(), " associations after filtering invalid values")
#   # Pick the lead variant from GWAS data
#   lead_var_gwas = dplyr::arrange(gwas_data, p_value)[1,]
# 
#   # Pick significant associations in eQTL Catalogue with lead GWAS variant ID for RNA-seq and Microarray studies and merge two tables
#   # "ROSMAP"
#   # "Schmiedel_2018",
#   rnaseq_study_names <- c("Alasoo_2018", "BLUEPRINT", "BrainSeq", "GENCORD", "GEUVADIS",
#                           "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016",
#                           "Schwartzentruber_2018", "TwinsUK", "van_de_Bunt_2015")
#   rnaseq_studies_sign_assocs <- get_significant_assocs_of_var(study_ids = rnaseq_study_names,
#                                                               variant_id = lead_var_gwas$variant_id)
# 
#   microarray_study_names <- c("CEDAR", "Fairfax_2014", "Kasela_2017", "Naranbhai_2015", "Fairfax_2012")
#   microarray_studies_sign_assocs <- get_significant_assocs_of_var(study_ids = microarray_study_names,
#                                                                   variant_id = lead_var_gwas$variant_id,
#                                                                   quant_method = "microarray")
# 
#   all_studies_sign_assocs <- rnaseq_studies_sign_assocs %>% rbind(microarray_studies_sign_assocs)
# 
# 
#   ## Plot pvalues of each significant association
#   # fetch all significant transcript QTLs
#   sign_tx_rs4239702 <-  get_significant_assocs_of_var(study_ids = rnaseq_study_names,
#                                                       variant_id = lead_var_gwas$variant_id,
#                                                       quant_method = "tx")
# 
#   all_studies_sign_assocs_with_tx <- all_studies_sign_assocs %>% rbind(sign_tx_rs4239702)
#   sign_assocs_to_plot_scatter <- all_studies_sign_assocs_with_tx %>% group_by(study_id, qtl_group) %>%
#     dplyr::filter(pvalue == min(pvalue)) %>%
#     summarise(condition_label = condition_label,
#               molecular_trait_id=molecular_trait_id,
#               pvalue = pvalue,
#               tissue_label = tissue_label,
#               quant_method = quant_method) %>%
#     arrange(pvalue) %>%
#     mutate(context = paste0(study_id,"_",condition_label))
# 
#   sign_assocs_to_plot_scatter$quant_method[sign_assocs_to_plot_scatter$quant_method=="tx"] <- "RNA-seq Transcript usage"
#   sign_assocs_to_plot_scatter$quant_method[sign_assocs_to_plot_scatter$quant_method=="ge"] <- "RNA-seq Gene expression"
#   sign_assocs_to_plot_scatter$quant_method[sign_assocs_to_plot_scatter$quant_method=="microarray"] <- "Microarray Gene expression"
#   sign_assocs_to_plot_scatter$context <- factor(sign_assocs_to_plot_scatter$context, levels = rev(unique(sign_assocs_to_plot_scatter$context)))
# 
#   sign_assocs_to_plot_scatter$tissue_label <- factor(sign_assocs_to_plot_scatter$tissue_label, levels = unique(sign_assocs_to_plot_scatter$tissue_label))
# 
#   base_plot <- ggplot(sign_assocs_to_plot_scatter, aes(x=-log10(pvalue), y = context, color = tissue_label, shape = quant_method))
#   final_plot <- base_plot + geom_point() + theme_bw()+ scale_color_brewer(palette="Dark2") +
#     xlab(expression(paste("-",log[10], " p-value"))) +
#     ylab("Study and Condition") +
#     labs(color = "Cell Type", shape = "Quantification method")
# 
#   save_ggplots(plot = final_plot, path = "eQTL_API_usecase_figures", filename = "sign_eqtl_scatter_shape", height = 10, width = 18)
#   final_plot
# 
# 
#   # Get all associations in cis region of each significant QTL
#   all_coloc_data <- list()
#   all_coloc_data$eqtl_assocs_in_region <- apply(all_studies_sign_assocs, 1, get_assocs_in_region)
# 
# 
#   # Perform colocalisation analysis.
#   all_coloc_data$colocs <- sapply(all_coloc_data$eqtl_assocs_in_region, get_colocs, gwas_data_for_trait=gwas_data)
# 
# 
#   # Prepare merged facet plottable dataframe
#   plottable_merged_data <- make_assocs_list_to_merged_plottable(all_coloc_data)
# 
# 
#   # Plot faceted manhattan figure with multiple
#   all_plots_faceted <- plot_faceted_multi_manhattans(merged_eqtl_assocs = plottable_merged_data,
#                                                      gwas_data_for_trait = gwas_data,
#                                                      save_plot = TRUE,
#                                                      save_dir = "eQTL_API_usecase_figures",
#                                                      save_filename = "merged_manhattan", no_GWAS_plot = TRUE)
#   all_plots_faceted
# 
# 
#   # Plot specific figures in faceted plot
#   plottable_merged_data_filt <- plottable_merged_data %>%
#     dplyr::filter(study_id %in% c("BLUEPRINT", "Quach_2016", "CEDAR", "Fairfax_2014")) %>%
#     dplyr::filter(qtl_group %in% c("monocyte_naive", "monocyte", "monocyte_CD14")) %>%
#     dplyr::filter(molecular_trait_id %in% c("ENSG00000101017", "ILMN_1779257"))
# 
#   filt_plots_faceted <- plot_faceted_multi_manhattans(merged_eqtl_assocs = plottable_merged_data_filt,
#                                                       gwas_data_for_trait = gwas_data,
#                                                       save_plot = TRUE,
#                                                       save_dir = "eQTL_API_usecase_figures",
#                                                       save_filename = "merged_manhattan_filt", no_GWAS_plot = TRUE)
#   filt_plots_faceted
# 
# 
# 
#   ## Analyse transcript usage associations in Alasoo_2018
#   sign_tx_Alasoo_2018 <- sign_tx_rs4239702 %>% dplyr::filter(study_id=="Alasoo_2018")
# 
#   all_coloc_data_tx <- list()
#   all_coloc_data_tx$eqtl_assocs_in_region <- apply(sign_tx_Alasoo_2018, 1, get_assocs_in_region)
# 
#   all_coloc_data_tx$colocs <- sapply(all_coloc_data_tx$eqtl_assocs_in_region, get_colocs, gwas_data_for_trait=gwas_data)
#   merged_plottable_df <-  make_assocs_list_to_merged_plottable(all_coloc_data_tx, quant_method = "tx")
# 
#   manhattan_tx_gwas <- plot_faceted_multi_manhattans(merged_eqtl_assocs = merged_plottable_df,
#                                                      gwas_data_for_trait = gwas_data,
#                                                      save_plot = TRUE,
#                                                      save_dir = "eQTL_API_usecase_figures",
#                                                      save_filename = "merged_manhattan_tx")
#   manhattan_tx_gwas
# 
# 
#   ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=78)
# 
#   CD_40_ENST00000466205_exons <- getBM(attributes=c('ensembl_exon_id','exon_chrom_start','exon_chrom_end',
#                                                     'strand', 'chromosome_name'),
#                                        filters = 'ensembl_transcript_id',
#                                        values ="ENST00000466205",
#                                        mart = ensembl)
# 
#   CD_40_ENST00000466205_exons$strand <- "+"
# 
#   exons_grange <- makeGRangesFromDataFrame(CD_40_ENST00000466205_exons, ignore.strand = FALSE, seqnames.field = "chromosome_name", start.field = 'exon_chrom_start', end.field = 'exon_chrom_end', strand.field = 'strand', keep.extra.columns = TRUE)
#   CD_40_exons_list <- GRangesList()
#   CD_40_exons_list[['ENST00000466205']] <- exons_grange
# 
#   region_coords = c(45980000, 46200000)
#   CD40_exons_plot <- plotTranscripts(exons = CD_40_exons_list, rescale_introns = FALSE, region_coords = region_coords)
#   CD40_exons_plot
# 
# 
#   joint_plot = cowplot::plot_grid(manhattan_tx_gwas, filt_plots_faceted, CD40_exons_plot,
#                                   align = "v", ncol = 1, rel_heights = c(2,4, 1))
# 
#   save_ggplots(plot = joint_plot, path = "eQTL_API_usecase_figures", filename = "new_tx_merged_plot_with_exon", height = 28, width = 15)
#   joint_plot
# }

