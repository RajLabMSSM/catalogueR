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

  



max_file_size <- function(top_dir, recursive=T){
  files <- list.files(top_dir, full.names = T, recursive = recursive)
  inf <- file.info(files) 
  print(paste("File size info (Mb) for",length(files),"files."))
  summary(inf$size/1e6)
}


timeit <- function(func, digits=3){
  start = Sys.time()
  out <- func
  end = Sys.time()
  print(paste("Completed in",round(end-start,digits = digits),"seconds."))
  return(out)
}

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
                             bp_lower=NULL,
                             multithread_tabix=F){
  if(use_tabix){
    # Tabix is about ~17x faster than the REST API.
    gwas.qtl <- catalogueR.fetch_tabix(unique_id=unique_id,
                                       quant_method="ge",
                                       infer_region=T,
                                       gwas_data=gwas_data,
                                       is_gwas=F,
                                       nThread = ifelse(multithread_tabix,nThread,1))
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
                                                    multithread_tabix=F,
                                                    nThread=4,
                                                    split_files=T,
                                                    merge_with_gwas=F,
                                                    force_new_subset=F,
                                                    progress_bar=T,
                                                    genome_build="hg19"){
  # WARNING!: pbmclapply only worked on Linux/Unix systems (e.g. mac) and NOT on Windows.
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
    split_path <- file.path(output_path, qtl_id, paste0(loc,"_locus___",qtl_id,".tsv.gz"))
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
                                       nThread=ifelse(multithread_tabix,nThread,1),
                                       multithread_tabix=multithread_tabix,
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
      if(split_files){data.table::fwrite(gwas.qtl, split_path, sep="\t", 
                                         nThread = ifelse(multithread_tabix,nThread,1))}
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
                           multithread_tabix=F,
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
                                                             genome_build=genome_build,
                                                             multithread_tabix=multithread_tabix)
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



catalogueR.gather_top_eVariants <- function(root_dir="/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019", 
                                            save_path="./eQTL_catalogue_topHits.tsv.gz",
                                            nThread=4,
                                            criterion="top_eVariant",
                                            pval_thresh=1e-5){
  library(dplyr)
  # root_dir = "/Volumes/Steelix/eQTL_catalogue/Nalls23andMe_2019"; nThread=4;
  qtl.paths <- list.files(root_dir, 
                          pattern="*.tsv.gz$", full.names = T, recursive = T) 
  # Group QTLs so you can get some sense of progress
  qtl.groups <- unique(basename(dirname(qtl.paths)))
  
  topQTL <- timeit(
    lapply(qtl.groups, function(group){
      message("+ ",group)
      qtl.paths_select <- qtl.paths[basename(dirname(qtl.paths))==group]
      
      topqtl <- parallel::mclapply(qtl.paths_select, function(x){
        message(basename(x))
        top_qtls <- NULL
        try({
          qtl <- data.table::fread(x, nThread = 1)
          print(dim(qtl))
          qtl.id <- strsplit(gsub(".tsv.gz$","",basename(x)), split = "___")[[1]][2]
          qtl <- qtl %>% 
            dplyr::mutate(eGene = ifelse((!is.na(gene.QTL) & gene.QTL!=""),gene.QTL,molecular_trait_id.QTL),
                          qtl.ID = qtl.id)  
         
            if(criterion=="pval_thresh"){
              # Get sig eVariants
              top_qtls <- subset(qtl, pvalue.QTL < pval_thresh)
            } else if(criterion=="top_eVariant"){
              # Get top eVariants per eGene
              top_qtls <- qtl %>% 
                dplyr::group_by(eGene) %>%  
                dplyr::top_n(n = 1, wt = -pvalue.QTL) %>% 
                data.table::data.table()
            }
          
          
        }) 
        return(top_qtls)
      }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)
      return(topqtl)
    }) %>% data.table::rbindlist(fill=T)
    
  ) 
  
  if(save_path!=F){
    dir.create(dirname(save_path),showWarnings = F, recursive = T)
    data.table::fwrite(topQTL, save_path, nThread = nThread, sep="\t")
  } 
  return(topQTL)
}


catalogueR.top_eVariants_overlap <- function(topQTL, 
                                             gwas_dataset="./Data/GWAS/Nalls23andMe_2019",
                                             gwas_min_support=1,
                                             qtl_pvalue_thresh=NULL){
  # topQTL <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_catalogue_topHits.tsv.gz", nThread = 4)
  merged_DT <- merge_finemapping_results(dataset = gwas_dataset,
                                         minimum_support = gwas_min_support,
                                         include_leadSNPs = T,
                                         xlsx_path = F)
  merged_DT$Locus <- merged_DT$Gene
  qtl.cols <- grep(".QTL$",colnames(topQTL), value = T)
  topQTL <- data.table:::merge.data.table(x=merged_DT,
                                          y=subset(topQTL, select=c("Locus",'SNP',"qtl.ID", qtl.cols, "eGene")),
                                          by=c("Locus","SNP"),
                                          all=T)
  
  if(is.null(qtl_pvalue_thresh)){
    qtl_pvalue_thresh <- 5e-8 / length(unique(topQTL$qtl.ID)) /length(unique(topQTL$eGene))
  }
  
  qtl_sig <- topQTL %>%
    subset(!is.na(pvalue.QTL) & pvalue.QTL<p_thresh) %>%
    dplyr::group_by(qtl.ID, Locus) %>% top_n(n=1, wt = -pvalue.QTL) %>% 
    dplyr::mutate(Locus.eGene = paste0(Locus,"  (",eGene,")")) %>%
    data.table::data.table()
  return(qtl_sig)
  
}


catalogueR.file_names <- function(topQTL, root_dir="./"){
  printer("Constructing file names for each GWAS-QTL locus...") 
  # Find which files are a good place to start with coloc
  # select_loci <- unique(subset(topQTL, Consensus_SNP | Support>0 | leadSNP, select=c("Locus","qtl.ID"))$Locus)
  # subset(topQTL, Locus %in% select_loci)
  if("qtl.id" %in% colnames(topQTL)){topQTL$qtl.ID <- topQTL$qtl.id }
  if("Locus.GWAS" %in% colnames(topQTL)){topQTL$Locus <- topQTL$Locus.GWAS }
  top_files <- ( topQTL%>%
                  dplyr::mutate(file= file.path(qtl.ID,paste0(Locus,"_locus___",qtl.ID,".tsv.gz")) ))$file %>% unique()
  printer("+",length(top_files),"file names returned.")
  return(file.path(root_dir,top_files))
}


catalogueR.plot_top_eVariants_overlap <- function(topQTL,
                                                 save_path="./topQTL_allLoci.png",
                                                 qtl_pvalue_thresh=NULL){ 
  # topQTL <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_catalogue_topHits.tsv.gz", nThread = 4)
  if(is.null(qtl_pvalue_thresh)){
    qtl_pvalue_thresh <- 5e-8 / length(unique(topQTL$qtl.ID)) /length(unique(topQTL$eGene))
  }
  topQTL_sig <- subset(topQTL, pvalue.QTL<qtl_pvalue_thresh)
  
  max_x <- length(unique(topQTL_sig$Locus.eGene))
  
  locus.gene_plot <- ggplot(topQTL_sig, aes(x=Locus.eGene, y=qtl.ID)) + 
    # geom_raster(aes(fill=topQTL.prop)) +  
    geom_raster(aes(fill=pvalue.QTL)) + 
    scale_fill_gradient(low = "red", high = "blue", na.value = "transparent") +
    # Consensus SNP markers
    geom_point(data = subset(topQTL_sig, Consensus_SNP), aes(x=Locus.eGene, y=qtl.ID), color="goldenrod2", show.legend = T) +
    
    # UCS SNP markers
    geom_point(data = subset(topQTL_sig, Support>0), aes(x=Locus.eGene, y=qtl.ID), size=3, shape=1, color="green2", show.legend = T) +
    
    # lead GWAS SNP markers
    geom_point(data = subset(topQTL_sig, leadSNP), aes(x=Locus.eGene, y=qtl.ID), size=5, shape=1, color="black", show.legend=T) +
    coord_cartesian(xlim=c(0,max_x+2), clip = 'off') +
    annotate("text", x=max_x+2, y=3, label = "• Consensus SNP", vjust = 2, hjust = 0, colour = "goldenrod2", size=3) +
    annotate("text",x=max_x+2, y=3, label = "○ Union Credible Set", vjust = 0, hjust = 0, colour = "green2", size=3) +
    annotate("text", x=max_x+2, y=3, label = "○ Lead GWAS SNP", vjust = -2,  hjust = 0, colour = "black", size=3) +
    scale_y_discrete(drop = F) +
    # facet_grid(facets = .~ Locus.eGene, scales = "free_x") + 
    scale_x_discrete(position = "top") +
    theme_bw() +
    labs(title=paste("Top eVariants per QTL dataset per GWAS Locus\n",
                     "Overlapping lead GWAS SNPs only"),
         # subtitle="Overlapping lead GWAS SNPs only",
         x="\nGWAS Locus (QTL eGene)\n") +
    theme(plot.title = element_text(hjust = .5), 
          plot.subtitle = element_text(hjust = .5),
          axis.text.x = element_text(angle=45, hjust=0))
  print(locus.gene_plot)
  
  if(save_path!=F){
    # save_path <- "./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/topQTL_leadGWAS-SNPs.png"
    ggsave(save_path, plot=locus.gene_plot, height=10, width=10)
  } 
  return(locus.gene_plot)
}



# #### Check QTL overlap
# 
# ```{r}
# gwas.qtl <- catalogueR.merge_gwas_qtl(gwas_data = merged_DT,
#                                       qtl.subset = qtl.dat) 
# qtl.sig <-  subset(gwas.qtl, pvalue.QTL<.05 & (Consensus_SNP|leadSNP) ) 
# 
# 
# ggplot(data=qtl.sig, aes(x=SNP, y=qtl.id, fill=-log10(pvalue.QTL))) +   
#   geom_raster() + 
#   # geom_point(data=subset(top.snps_melt, leadQTL), shape="*", color="cyan", size=5) +
#   facet_grid(facets = .~gene.QTL + Locus, scales = "free_x") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_continuous(low="thistle2", high="darkred", 
#                         guide="colorbar",na.value="transparent") 
# 
# # QTL pvalue distribution plots
# ggplot(data=gwas.qtl, aes(x=POS, y=-log10(pvalue.QTL), color=-log10(pvalue.QTL))) + 
#   geom_point(size=.1) +
#   geom_hline(yintercept = -log10(5e-8), show.legend = F, alpha=.5, size=.1) +
#   facet_grid(facets = Locus + qtl.id ~ gene.QTL, 
#              scales = "free") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         strip.text.y = element_text(angle = 0, margin = margin(0, 0, 0, 0, "cm")))  
# ```


###############----------COLOCALIZATION---------------##################

COLOC.report_summary <- function(coloc.res, PP_threshold=.8){ 
  # MAF = dataset1$MAF) 
  hypothesis_key <- setNames(
    c("Neither trait has a genetic association in the region.",
      "Only trait 1 has a genetic association in the region.",
      "Only trait 2 has a genetic association in the region.",
      "Both traits are associated, but with different causal variants (one in each dataset).",
      "Both traits are associated and share a single causal variant.") ,
    c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")) 
  # Report hypothess results
  printer("Hypothesis Results @ PP_threshold =",PP_threshold,":")
  true_hyp <-""
  for(h in names(hypothesis_key)){
    if(!is.na(coloc.res$summary[h])){
      if(coloc.res$summary[h]>=PP_threshold){
        hyp <- hypothesis_key[h]
        printer("    ",h,"== TRUE: **",hyp )
        true_hyp <- paste0(names(hyp),": ", hyp)
      }
    } else{
      printer("    ",h,"== FALSE: ")
    } 
  } 
  
  # Save raw results   
  coloc_DT <- coloc.res$results
  # Process results  
  coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
  colocalized_snps <- subset(coloc_DT, Colocalized==T)$snp# subset(coloc_DT, Colocalized==1)$SNP
  subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps,sep=", "))
  if(!is.na(coloc.res$summary)["PP.H4.abf"] ){
    if((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) & 
       (coloc.res$summary["PP.H4.abf"]/coloc.res$summary["PP.H3.abf"] >= 2)){
      # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
      report <- paste("Datasets colocalized")  
    } else {report <- paste("Datasets NOT colocalized") }
  } else { report <- paste("Datasets NOT colocalized")}   
  printer("+ COLOC::",report,"at: PP.H3 + PP.H4 >=",PP_threshold," and PP.H3 / PP.H4 >= 2.") 
  return(coloc_DT)
}



catalogueR.get_colocs <- function(qtl.egene, 
                                  gwas.region,
                                  merge_by_rsid=T,
                                  PP_threshold=.8,
                                  verbose=T){
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
    if(verbose){  message("catalogueR:COLOC:: No SNPs shared between GWAS and QTL subsets.")  }
  
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
    
    # wrap <- ifelse(verbose, function(x)x, suppressMessages)
    coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, 
                                 dataset2 = gwas_dataset,
                                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) # defaults
    coloc_res$Locus <- gwas_shared$Locus[1]
    if(verbose){
      report <- COLOC.report_summary(coloc.res = coloc_res, 
                                     PP_threshold = PP_threshold)
    }
  
    }
  return(coloc_res)
}




catalogueR.run_coloc <- function(gwas.qtl_paths,
                                 save_path="./coloc_results.tsv.gz",
                                 nThread=3,
                                 top_snp_only=T,
                                 split_by_group=F){
  # gwas.qtl_paths <- list.files("/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019", recursive = T, full.names = T)
  # gwas.qtl_paths <- list.files("../eQTL_catalogue/Nalls23andMe_2019", recursive = T, full.names = T)
  # gwas.qtl_paths <- list.files("/Volumes/Steelix/eQTL_catalogue/Nalls23andMe_2019", recursive = T, full.names = T)
  # gwas.qtl_paths <- file.path("/Volumes/Steelix/eQTL_catalogue/Nalls23andMe_2019",top_files)
  # qtl.ID=unique(qtl.dat$qtl.id)[1]; eGene=unique(qtl.dataset$gene.QTL)[1];
  # qtl.ID = "Fairfax_2014.monocyte_naive"; eGene = "TSC22D2"; qtl.path=gwas.qtl_paths[1]
  # qtl.path="/Volumes/Scizor/eQTL_catalogue/Nalls23andMe_2019/Alasoo_2018.macrophage_IFNg/MCCC1_locus__&__Alasoo_2018.macrophage_IFNg.tsv.gz"
  qtl.groups <- unique(basename(dirname(gwas.qtl_paths)))
  
  # Iterate over QTL groups 
  coloc_QTLs <- timeit(lapply(qtl.groups, function(group){
    message("+ QTL Group = ",group)
    gwas.qtl_paths_select <- gwas.qtl_paths[basename(dirname(gwas.qtl_paths))==group]
    
    # ---- Iterate over QTL datasets 
    coloc_qtls <- parallel::mclapply(gwas.qtl_paths_select, function(qtl.path){
      qtl.ID <- strsplit(gsub(".tsv.gz","", basename(qtl.path)), "___")[[1]][2]
      gwas.locus <- strsplit(gsub(".tsv.gz","", basename(qtl.path)), "___")[[1]][1] 
      coloc_eGenes <- data.table::data.table()
      try({
        qtl.dat <- data.table::fread(qtl.path)
        message("++ GWAS =",gwas.locus," x ",length(unique(qtl.dat$gene.QTL))," eGenes")
        if(!"qtl.id" %in% colnames(qtl.dat)){qtl.dat <- cbind(qtl.dat, qtl.id=qtl.ID)}
        qtl.dataset <- subset(qtl.dat, qtl.id==qtl.ID & !is.na(gene.QTL) & gene.QTL!="")
        remove(qtl.dat)
        if("Effect" %in% colnames(qtl.dataset)){
          gwas_cols <- c("Locus","Locus.GWAS", "SNP", "CHR","POS", "P", "Effect", "StdErr", "Freq", "MAF", "N_cases", "N_controls", "proportion_cases", "A1", "A2")
          gwas_cols <- gwas_cols[gwas_cols %in% colnames(qtl.dataset)]
          gwas.region <- subset(qtl.dataset, select=gwas_cols)
        }  
        
        # ---- Iterate over QTL eGenes
          coloc_eGenes <- parallel::mclapply(unique(qtl.dataset$gene.QTL), function(eGene){ 
            printer("+++ QTL eGene =",eGene)
            qtl.egene <- subset(qtl.dataset, gene.QTL==eGene) 
            # gwas.region <- subset(gwas_data, SNP %in% qtl.egene$rsid.QTL) 
            coloc_res <- catalogueR.get_colocs(qtl.egene = qtl.egene, 
                                               gwas.region = gwas.region, 
                                               merge_by_rsid = T,
                                               verbose = F)
            coloc_summary <- as.list(coloc_res$summary)
            coloc_results <- coloc_res$results
            if(top_snp_only){
              coloc_results <- (coloc_results %>% top_n(n=1, wt=SNP.PP.H4))[1,]
            }
            
            
            coloc_DT <- data.table::data.table(Locus.GWAS=coloc_res$Locus, 
                                               qtl.id=qtl.ID,
                                               eGene=eGene,
                                               coloc_results,
                                               PP.H0=coloc_summary$PP.H0.abf,
                                               PP.H1=coloc_summary$PP.H1.abf,
                                               PP.H2=coloc_summary$PP.H2.abf,
                                               PP.H3=coloc_summary$PP.H3.abf,
                                               PP.H4=coloc_summary$PP.H4.abf)
            return(coloc_DT)
          }) %>% data.table::rbindlist(fill=T)   
        }) # end try()
        return(coloc_eGenes)  
      }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)  
    
    if(split_by_group){
      split_path <- file.path(dirname(save_path),group)
      printer("Saving split file ==>", split_path) 
      dir.create(dirname(split_path), showWarnings = F, recursive = T)
      data.table::fwrite(coloc_qtls, file = split_path)
      return(split_path)
    } else { return(coloc_qtls)}
   
  }) ) # end timeit()
  
  if(split_by_group){
    printer("+ Returning split file paths.")
  } else {
    coloc_QTLs <- data.table::rbindlist(coloc_QTLs, fill=T)
    # Save all coloc results in one dataframe
    if(save_path!=F){
      # save_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/coloc_topQTL.tsv"
      dir.create(dirname(save_path), showWarnings = F, recursive = T)
      data.table::fwrite(coloc_QTLs, file = save_path, sep='\t')
    } 
  } 
  return(coloc_QTLs) 
}




catalogueR.plot_coloc_summary <- function(coloc_QTLs, 
                                          topQTL,
                                          PP_thresh = .8, 
                                          save_path="./coloc_topQTL.png"){ 
  # topQTL <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_catalogue_topHits.tsv.gz", nThread = 4)
  # coloc_QTLs <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz", nThread = 4)
  # PP_thresh <- .99
 
  coloc_plot <- subset(coloc_QTLs, !is.na(Locus.GWAS)) %>% 
    dplyr::mutate(PP.H4.thresh = ifelse(PP.H4>=PP_thresh, PP.H4,NA),
                  PP.Hyp4= ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4/PP.H3 >= 2), PP.H4,NA)) %>%
    tidyr::separate(col = qtl.id, into = c("QTL.group","id"), sep = "\\.", remove = F) %>%
    # only plot colocalized loci-egene combinations (otherwise, plot will be massive)
    data.table::data.table() %>%
    data.table:::merge.data.table(y=topQTL,
                                  by.x=c("Locus.GWAS","eGene","qtl.id"),
                                  by.y=c("Locus","eGene","qtl.ID"), 
                                  all = T)  %>%
    # subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh) | (Support>0) | (leadSNP)), .drop=F) %>%
    subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop=F) %>%
    dplyr::mutate(Locus.eGene = paste0(Locus.GWAS," (",eGene,")"),
                  Consensus_SNP = Consensus_SNP,
                  Union_Credible_Set = Support>0,
                  Lead_GWAS_SNP = leadSNP) %>%
    # Get only the top eGene per locus
    dplyr::group_by(Locus.GWAS, qtl.id) %>% 
    top_n(n=1,wt=-pvalue.QTL) %>% 
    data.table::data.table()
  
  max_x <- length(unique(coloc_plot$Locus.eGene))
   
  # raster plot 
  gg_coloc <- ggplot(data=coloc_plot, aes(x=Locus.eGene, y=qtl.id, fill=PP.Hyp4)) + 
    geom_raster() +   
    # scale_fill_continuous(limits=c(minPP4,1)) +
    scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue",.7), high = scales::alpha("red",.7)) +
    scale_y_discrete(drop=F) +
    scale_x_discrete(position = "top") +
    # scale_fill_manual(limits=c(0,1), palette="Spectral") +
    # facet_grid(facets = . ~ Locus.GWAS + eGene, 
    #            switch = "y",scales = "free") +
    labs(title=paste0("Colocalized GWAS loci x QTL loci (>",PP_thresh*100,"% probability)\n"),
         x="GWAS Locus (QTL eGene)",
         fill="Colocalization\nProbability") +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5),
          axis.text.x = element_text(angle=45, hjust=0), 
          panel.border = element_rect(colour = "transparent")) 
    # scale_color_manual(name = "SNP group x lead eVariant overlap",
    #               breaks = c("Consensus_SNP", "Union_Credible_Set", "Lead_GWAS_SNP"),
    #               values = c(Consensus_SNP = "goldenrod2", Union_Credible_Set = "green2", Union_Credible_Set = "grey30") ) 
    # coord_cartesian(xlim=c(0,max_x+2), clip = 'off') +
    # annotate("text", x=max_x+2, y=3, label = "• Consensus SNP", vjust = 2, hjust = 0, colour = "goldenrod2", size=3) +
    # annotate("text",x=max_x+2, y=3, label = "○ Union Credible Set", vjust = 0, hjust = 0, colour = "green2", size=3) +
    # annotate("text", x=max_x+2, y=3, label = "○ Lead GWAS SNP", vjust = -2,  hjust = 0, colour = "black", size=3)
  # Consensus SNP markers
  # geom_point(data = subset(coloc_plot, Consensus_SNP), aes(x=Locus.eGene, y=qtl.id, color="Consensus_SNP"), color="goldenrod2", show.legend = T) +
  #   # UCS SNP markers
  #   geom_point(data = subset(coloc_plot, Union_Credible_Set), aes(x=Locus.eGene, y=qtl.id, color="Union_Credible_Set"), size=3, shape=1, color="green2", show.legend = T) +
  #   # lead GWAS SNP markers
  #   geom_point(data = subset(coloc_plot, Lead_GWAS_SNP), aes(x=Locus.eGene, y=qtl.id, color="Lead_GWAS_SNP"), size=5, shape=1, color="black", show.legend=T) +
  print(gg_coloc)
  
  if(save_path!=F){
    # save_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/coloc_PP80.png"
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    ggsave(save_path,
           plot=gg_coloc, height = 11, width=25)
  }
  # Return merged data
 return(coloc_plot)
}










# Insanely large faceted manhattan plots of QTL datasets
catalogueR.manhattan_facets <- function(gwas.qtl_paths,
                                        qtl.thresh = 1e-5,
                                        gwas_label="GWAS"){
  # gwas.qtl_paths <- list.files("/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019", full.names = T, recursive = T)
  # gwas.qtl_paths <- gwas.qtl_paths[1:10]
  PP.H4.thresh=.8
  coloc_QTLs_sig <- coloc_QTLs %>% dplyr::mutate(PP.H4.thresh = ifelse(PP.H4>=PP_thresh, PP.H4,NA),
                PP.Hyp4= ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4/PP.H3 >= 2), PP.H4,NA)) %>% 
    subset(!is.na(PP.Hyp4))# %>% 
    # subset(Locus.GWAS %in% c("BIN3","MED12L","LRRK2"))
  top_files <- catalogueR.file_names(topQTL = coloc_QTLs_sig,
                                    root_dir = "/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019") 
  qtl.dat <- rbind.file.list(top_files)
  
  qtl.dat$qtl.id <- gsub(".tsv.gz$","",strsplit(basename(gwas.qtl_paths),"___")[[1]][2] )
  if(!"Support" %in% colnames(qtl.dat)){
    qtl.dat <- find_consensus_SNPs(qtl.dat, consensus_thresh = 2)
  }
  plot.subset <- qtl.dat %>% dplyr::mutate(UCS = Support>0, 
                                        GWAS_label=gwas_label,
                                        MB=POS/1000000) %>% 
    subset(pvalue.QTL<qtl.thresh) %>%
    add_eGene_col()  
    # Only show colocalized plots within EACH LOCUS.
    ## Otherwise, there's too many to see.
  
    gg.gwas <- ggplot(plot.subset, aes(x=MB, y=-log10(P), color=-log10(P))) + 
      geom_hline(yintercept = -log10(5e-8), alpha=.5, linetype = "dashed", size=.5) + 
      geom_point(size=.5) + 
      labs(x=NULL) +
      facet_grid(facets = GWAS_label ~ Locus , 
                 scales = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_blank(), 
            strip.text.y = element_text(angle = 0)) 
    
    gg.qtls <-  ggplot(plot.subset, aes(x=MB, y=-log10(pvalue.QTL))) +  
      geom_point(aes(color=eGene), size=.5, show.legend = F) + 
      geom_hline(yintercept = -log10(qtl.thresh), alpha=.5, linetype = "dashed", size=.5) + 
      facet_grid(facets = qtl.id ~ eGene + Locus, 
                 scales = "free_x") +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
    
    # Merge plots
  gg.gwas + gg.qtls + 
      patchwork::plot_layout(ncol = 1)#, heights = c(.1,1))
  
  
  
  
  
  # Group and melt
  gwas.qtl.melt <-  setDT(gwas.qtl)[, .(Count = uniqueN(SNP[UCS==T],na.rm = T)),
                                    by=c("Locus","qtl.id","SNP","gene.QTL")] #"leadSNP","Consensus_SNP","UCS"
  gwas.qtl.melt <- gwas.qtl.melt %>% 
    dplyr::group_by(Locus, qtl.id, SNP) %>% 
    dplyr::summarise(Count = sum(Count, na.rm=T),
                     Genes = paste(gene.QTL,collapse=", "))
  
  gwas.qtl.melt[gwas.qtl.melt$Count==0,"Count"] <- NA
  gg_gwas.qtl <- ggplot(data=gwas.qtl.melt, aes(x=qtl.id, y=Locus, fill=Count)) +
    geom_raster() +  
    # scale_fill_manual(values = consensus_colors) +  
    scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue",.7), high = scales::alpha("red",.7)) +
    geom_point(aes(size=ifelse(Count>0, "dot", "no_dot")), show.legend = F, alpha=.8, color="white") +
    scale_size_manual(values=c(dot=.5, no_dot=NA), guide="none") +
    labs(fill = "UCS SNP Count") +
    theme_classic() +
    theme(legend.position = "top",  
          legend.title.align = .5,
          axis.text.x = element_text(angle = 40, hjust = 1),
          # legend.background =  element_rect(fill = "lightgray"),
          legend.key = element_rect(colour = "gray60"), 
          legend.text = element_text(size = 8),
          legend.text.align = .5,
          # legend.key.size = unit(.5, units = "cm" ),
          legend.box="horizontal",
          # panel.background = element_rect(fill = 'transparent'),
          # panel.grid = element_line(color="gray", size=5),
          panel.grid.major = element_line(color="grey", size=5) ) +  
    guides(color = guide_legend(nrow = 1, reverse = F,
                                title.position = "top",
                                # label.position = "top",
                                title.hjust = .5,
                                label.hjust = -1)) +
    # Keep unused levels/Loci
    scale_y_discrete(drop=FALSE)
  print(gg_gwas.qtl)
  return(gg_gwas.qtl)
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

find_consensus_SNPs <- function(finemap_DT,
                                verbose=T,
                                credset_thresh=.95,
                                consensus_thresh=2,
                                sort_by_support=T,
                                exclude_methods=NULL){
  printer("+ Identifying Consensus SNPs...",v=verbose)
  exclude_methods <- append(exclude_methods,"mean")
  # Find SNPs that are in the credible set for all fine-mapping tools
  CS_cols <- colnames(finemap_DT)[endsWith(colnames(finemap_DT),".Credible_Set")]
  CS_cols <- CS_cols[!(CS_cols %in% paste0(exclude_methods,".Credible_Set"))]
  if(consensus_thresh=="all"){consensus_thresh<-length(CS_cols)}
  printer("++ support_thresh =",consensus_thresh)
  # Get the number of tools supporting each SNP
  ## Make sure each CS is set to 1
  support_sub <- subset(finemap_DT, select = CS_cols) %>% data.frame()
  support_sub[sapply(support_sub, function(e){e>1})] <- 1
  finemap_DT$Support <- rowSums(support_sub, na.rm = T)
  finemap_DT$Consensus_SNP <- finemap_DT$Support >= consensus_thresh
  # Sort
  if(sort_by_support){
    finemap_DT <- finemap_DT %>% arrange(desc(Consensus_SNP), desc(Support))
  }
  
  # Calculate mean PP
  printer("+ Calculating mean Posterior Probability (mean.PP)...")
  PP.cols <- grep(".PP",colnames(finemap_DT), value = T)
  PP.cols <- PP.cols[!(PP.cols %in% paste0(exclude_methods,".PP"))]
  PP.sub <- subset(finemap_DT, select=c("SNP",PP.cols)) %>% data.frame()# %>% unique()
  PP.sub[is.na(PP.sub)] <- 0
  if(NCOL(PP.sub[,-1]) > 1){
    finemap_DT$mean.PP <- rowMeans(PP.sub[,-1])
  } else{
    finemap_DT$mean.PP <- PP.sub[,-1]
  }
  finemap_DT$mean.Credible_Set <- ifelse(finemap_DT$mean.PP>=credset_thresh,1,0)
  
  # PP.sub %>% arrange(desc(mean.PP)) %>% head()
  printer("++",length(CS_cols),"fine-mapping methods used.")
  printer("++",dim(subset(finemap_DT,Support>0))[1],"Credible Set SNPs identified.")
  printer("++",dim(subset(finemap_DT,Consensus_SNP==T))[1],"Consensus SNPs identified.")
  return(finemap_DT)
}


add_eGene_col <- function(qtl.dat){
  # Remove any old col
  if("eGene" %in% colnames(qtl.dat)){qtl.dat <- subset(qtl.dat, select=-eGene)}
  qtl.dat <- qtl.dat %>% dplyr::mutate(eGene = ifelse(!is.na(gene.QTL) & gene.QTL!="",gene.QTL, molecular_trait_id.QTL))
  return(qtl.dat)
}
