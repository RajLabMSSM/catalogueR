
###############----------COLOCALIZATION---------------##################






#' Report coloc results
#' 
#' @family coloc
#' @keywords internal 
COLOC.report_summary <- function(coloc.res, 
                                 PP_threshold=.8){ 
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




#' Run coloc on GWAS-QTL object
#' 
#' @family coloc
#' @examples 
#' library(dplyr)
#' data("BST1__Alasoo_2018.macrophage_IFNg")
#' 
#' qtl.egene <- data.frame(BST1__Alasoo_2018.macrophage_IFNg)[,grep("*.QTL$|qtl_id|SNP",colnames(BST1__Alasoo_2018.macrophage_IFNg), value = T)]
#' sorted_egenes <- qtl.egene %>% dplyr::group_by(gene.QTL) %>% dplyr::summarise(mean.P=mean(pvalue.QTL), min.P=min(pvalue.QTL)) %>% dplyr::arrange(min.P)
#' qtl.egene <- subset(qtl.egene, gene.QTL==sorted_egenes$gene.QTL[1])
#' 
#' gwas.region <- data.frame(BST1__Alasoo_2018.macrophage_IFNg)[,grep("*.QTL$|qtl_id",colnames(BST1__Alasoo_2018.macrophage_IFNg), value = T, invert = T)]
#' coloc_res <- get_colocs(qtl.egene=qtl.egene, gwas.region=gwas.region)
get_colocs <- function(qtl.egene, 
                       gwas.region,
                       merge_by_rsid=T,
                       PP_threshold=.8,
                       method="abf",
                       verbose=T){
  # http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html
  # Subset to overlapping SNPs only
  
  # Check for duplicated SNPs
  qtl.egene <- qtl.egene[!duplicated(qtl.egene$SNP),]
  gwas.region <- gwas.region[!duplicated(gwas.region$SNP),]
  if(all(c("N_cases","N_controls") %in% colnames(gwas.region))){
    message("Inferring `N` (effective sample size) from `N_cases` and `N_controls`")
    gwas.region$N <- get_sample_size(gwas.region, effective_ss=T)$N
  }
  if(!"N" %in% colnames(gwas.region)){
    stop("`N` column (effective sample size) was not detected in gwas.region. Required for coloc analysis.")
  }
  
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
  
  #### Check MAF ####
  if(!"MAF" %in% colnames(gwas_shared)){
    warning("`MAF` column not provided in GWAS data. Borrowing MAF from QTL data instead.")
    gwas_shared$MAF <- eqtl_shared$maf.QTL
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
    eQTL_dataset = list(snp = eqtl_shared$variant_id,
                        pvalues = eqtl_shared$pvalue.QTL,
                        #If log_OR column is full of NAs then use beta column instead
                        beta = eqtl_shared$beta.QTL, 
                        N = (eqtl_shared$an.QTL)[1],#/2,  
                        MAF = eqtl_shared$maf.QTL, 
                        type = "quant"
                        )
    if("se.QTL" %in% colnames(eqtl_shared)){
      eQTL_dataset$varbeta <- eqtl_shared$se.QTL^2 
    }
    # GWAS data
    gwas_dataset = list(snp = gwas_shared$variant_id,
                        pvalues = gwas_shared$P,
                        #If log_OR column is full of NAs then use beta column instead
                        beta = gwas_shared$Effect,  
                        N = (gwas_shared$N)[1],#/2,   
                        MAF = gwas_shared$MAF,
                        type = "cc",  
                        s = 0.5  #This is actually not used, because we already specified varbeta above.
                        )
    if("StdErr" %in% colnames(gwas_shared)){
      gwas_dataset$varbeta <-gwas_shared$StdErr^2 
    }  
    
    # wrap <- ifelse(verbose, function(x)x, suppressMessages)
     
    # coloc::coloc.signals()
    if(method %in% c("single","cond","mask")){
      # The updated coloc requires varbeta in both datasets 
      color_res <- coloc::coloc.signals(dataset1 = eQTL_dataset, 
                                        dataset2 = gwas_dataset,
                                        method = "mask",
                                        p12 = 1e-5)
    }
    if(tolower(method)=="abf"){
      coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, 
                                   dataset2 = gwas_dataset,
                                   p12 = 1e-5) # defaults
    }
    
    coloc_res$Locus <- gwas_shared$Locus[1]
    if(verbose){
      report <- COLOC.report_summary(coloc.res = coloc_res, 
                                     PP_threshold = PP_threshold)
    }
    
  }
  return(coloc_res)
}





#' Iteratively run coloc on merged GWAS-QTL datatables
#' 
#' Runs colocalization tests (\code{\link{coloc::coloc.abf}}) on merged GWAS-QTL datatables
#' generated by \code{\link{catalogueR::eQTL_Catalogue.query}}.
#' 
#' Iterately runs coloc across each:
#' \itemize{
#' \item{QTL dataset}
#' \item{GWAS locus}
#' \item{QTL gene}
#' }
#' @return 
#' If \code{top_snp_only=T}, returns SNP-level stats for only the SNP 
#' with the highest colocalization probability (\emph{SNP.PP.H4})
#' If \code{top_snp_only=T}, returns SNP-level stats for every SNP.
#' In either case, summary-level coloc stats are added in the columns
#' \emph{PP.H0}, \emph{PP.H1}, \emph{PP.H2}, \emph{PP.H3}, \emph{PP.H4}.
#' @family coloc
#' @examples 
#' # With built-in data
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths() 
#' coloc_QTLs <- run_coloc(gwas.qtl_paths=gwas.qtl_paths, nThread=4, top_snp_only=T, save_path="~/Desktop/coloc_results.tsv.gz")
#' 
#' 
#' \dontrun{ 
#' # With full  Nalls et al data (not included)
#' gwas.qtl_paths <- list.files("/pd-omics/brian/eQTL_catalogue/Nalls23andMe_2019", recursive=T, full.names = T)
#' coloc_QTLs.Nalls2019 <- run_coloc(gwas.qtl_paths=gwas.qtl_paths[1:100], nThread=4, top_snp_only=T, save_path="~/Desktop/Nall2019.coloc_results.tsv.gz")
#' }
run_coloc <- function(gwas.qtl_paths,
                      save_path="./coloc_results.tsv.gz",
                      nThread=4,
                      top_snp_only=T,
                      split_by_group=F,
                      method="abf",
                      PP_threshold=.8,
                      gwas_sample_size=NULL){ 
  # gwas.qtl_paths <- list.files("~/Desktop/Mount_Sinai/COVID19_PheWAS/GWAS_queries", recursive = T, full.names = T) 
  # qtl.ID=unique(qtl.dat$qtl.id)[1]; eGene=unique(qtl.dataset$gene.QTL)[1];
  # qtl.ID = "Fairfax_2014.monocyte_naive"; eGene = "TSC22D2"; qtl.path=gwas.qtl_paths[1]
  # qtl.path="/Volumes/Scizor/eQTL_Catalogue/Nalls23andMe_2019/Alasoo_2018.macrophage_IFNg/MCCC1_locus__&__Alasoo_2018.macrophage_IFNg.tsv.gz"
  qtl.groups <- unique(unlist(lapply(gwas.qtl_paths, parse_gwas.qtl_path)))
  
  # Iterate over QTL groups 
  coloc_QTLs <- timeit(lapply(qtl.groups, function(group){
    message("+ QTL Group = ",group)
    gwas.qtl_paths_select <- gwas.qtl_paths[basename(dirname(gwas.qtl_paths))==group]
    
    # ---- Iterate over QTL datasets 
    coloc_qtls <- parallel::mclapply(gwas.qtl_paths_select, function(qtl.path,
                                                                     .gwas_sample_size=gwas_sample_size){
      qtl_ID <- parse_gwas.qtl_path(gwas.qtl_path = qtl.path, get_locus = F, get_qtl_id = T)
      gwas.locus <- parse_gwas.qtl_path(gwas.qtl_path = qtl.path, get_locus = T, get_qtl_id = F)
      coloc_eGenes <- data.table::data.table()
      try({
        qtl.dat <- data.table::fread(qtl.path) 
        message("++ GWAS = ",gwas.locus," x ",length(unique(qtl.dat$gene.QTL))," eGenes")
        #### Check column names ####
        
        ## Check sample size 
        if(!"N" %in% colnames(gwas.locus)){
          if(is.null(gwas_sample_size)){
            stop("`N` column (sample size) not detected. Please provide `sample_size=` argument instead.")
          }else{
            qtl.dat$N <- .gwas_sample_size
          } 
        } 
        if(!"qtl_id" %in% colnames(qtl.dat)) {qtl.dat <- cbind(qtl.dat, qtl_id=qtl_ID) }
        qtl.dataset <- subset(qtl.dat, qtl_id==qtl_ID & !is.na(gene.QTL) & gene.QTL!="")
        remove(qtl.dat) 
        
        # ---- Iterate over QTL eGenes
        coloc_eGenes <- lapply(unique(qtl.dataset$gene.QTL), function(eGene){ 
          printer("+++ QTL eGene =",eGene)
          # Extract QTL cols
          qtl.egene <- subset(qtl.dataset, gene.QTL==eGene)  
          # Extract the GWAS cols
          if("Effect" %in% colnames(qtl.egene)){
            gwas_cols <- c("Locus","Locus.GWAS", "SNP", "CHR","POS", "P", "Effect", "StdErr", "Freq",
                           "MAF","N","N_cases", "N_controls", "proportion_cases", "A1", "A2")
            gwas_cols <- gwas_cols[gwas_cols %in% colnames(qtl.egene)]
            gwas.region <- subset(qtl.egene, select=gwas_cols)
          }  
          # Run coloc
          coloc_res <- get_colocs(qtl.egene = qtl.egene, 
                                  gwas.region = gwas.region, 
                                  merge_by_rsid = T,
                                  PP_threshold = PP_threshold,
                                  method = method,
                                  verbose = F)
          coloc_summary <- as.list(coloc_res$summary)
          coloc_results <- coloc_res$results
          if(top_snp_only){
            coloc_results <- (coloc_results %>% top_n(n=1, wt=SNP.PP.H4))[1,]
          }
          # Rename columns
          colnames(coloc_results) <- gsub("\\.df1",".QTL",colnames(coloc_results))
          colnames(coloc_results) <- gsub("\\.df2",".GWAS",colnames(coloc_results))
          # Create results dataframe
          coloc_DT <- data.table::data.table(Locus.GWAS=coloc_res$Locus, 
                                             qtl_id=qtl_ID,
                                             gene.QTL=eGene,
                                             P=gwas.region$P,
                                             CHR=gwas.region$CHR,
                                             POS=gwas.region$POS,
                                             pvalue.QTL=qtl.egene$pvalue.QTL,
                                             coloc_results,
                                             PP.H0=coloc_summary$PP.H0.abf,
                                             PP.H1=coloc_summary$PP.H1.abf,
                                             PP.H2=coloc_summary$PP.H2.abf,
                                             PP.H3=coloc_summary$PP.H3.abf,
                                             PP.H4=coloc_summary$PP.H4.abf)
          return(coloc_DT)
        }) %>% data.table::rbindlist(fill=T)  # END ITERATE OVER eGenes
      }) # end try()
      return(coloc_eGenes)  
    }, mc.cores = nThread) %>% data.table::rbindlist(fill=T)  # END ITERATE OVER gwas.qtl_paths_select
    
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




# 
# add_old_coloc_results <- function(old_results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/old_results/COLOC_results_flip-gwasEffect.txt",
#                                   coloc_QTLs){
#   old_coloc <- data.table::fread(old_results_path)
#   old_coloc <- subset(old_coloc, !startsWith(Dataset2,"Fairfax")) %>% 
#     dplyr::rename(PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, 
#                   PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf,
#                   Locus.GWAS = Locus) 
#   old_coloc$eGene <- old_coloc$Locus.GWAS
#   old_coloc$qtl.id <- sub("^([^_]*_[^_]*)_", "\\1.",old_coloc$Dataset2)
#   old_coloc$qtl.id <- gsub("MESA_", "MESA.",old_coloc$qtl.id)
#   old_coloc$qtl.ID <- old_coloc$qtl.id
#   coloc_QTLs <- rbind(coloc_QTLs, old_coloc, fill=T)
#   return(coloc_QTLs)
# }




#' Plot summary of coloc tests
#' 
#' Plots the output from \code{catalogueR::run_coloc}.
#' @family coloc
#' @examples 
#' data("coloc_QTLs")
#' gg_coloc <- plot_coloc_summary(coloc_QTLs=coloc_QTLs, PP_thresh=.99, label_snp_groups=F)
plot_coloc_summary <- function(coloc_QTLs,   
                               PP_thresh=.8,  
                               label_snp_groups=F,
                               save_dir=".",
                               show_plot=T){    
  # merged_DT <- openxlsx::read.xlsx("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls23andMe_2019_results.xlsx")
  # merged_DT <- echolocatoR::update_CS_cols(merged_DT)   
  
  coloc_dat <- subset(coloc_QTLs, !is.na(Locus.GWAS)) %>%  
    dplyr::rename(qtl_ID=qtl_id) %>%
    dplyr::mutate(PP.H4.thresh = ifelse(PP.H4>=PP_thresh, PP.H4,NA),
                  PP.Hyp4 = ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4/PP.H3 >= 2), PP.H4,NA),
                  Locus.eGene = paste0(Locus.GWAS," (",gene.QTL,")")) %>%
    tidyr::separate(col = qtl_ID, into = c("group.QTL","id"), sep = "\\.", remove = F) %>%
    # only plot colocalized loci-egene combinations (otherwise, plot will be massive)
    subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop=F) %>%  
    data.table::data.table()
  if(nrow(coloc_dat)==0){
    stop("\n\n** No coloc tests reached the PPH4 threshold (",PP_thresh,").",
         "Try lowering `PP_thresh`. **")
  }
  
  needed_cols <- c("leadSNP","Consensus_SNP","Support")
  if(label_snp_groups){ 
    if(all(needed_cols %in% colnames(coloc_dat))){
      coloc_dat <- coloc_dat %>% 
        dplyr::group_by(Locus.GWAS,eGene,qtl_ID) %>% 
        dplyr::summarise(PP.Hyp4=max(PP.Hyp4),
                         leadGWAS.sigQTL = sum(leadSNP, na.rm=T),
                         Consensus.sigQTL = sum(Consensus_SNP, na.rm=T),
                         UCS.sigQTL = n_distinct(SNP[Support>0],na.rm = T)
        ) %>% 
        subset(!is.na(PP.Hyp4)) %>%
        data.table::data.table() 
    } else {
      missing_cols <- needed_cols[!needed_cols %in% colnames(coloc_dat)]
      warning("++ Cannot label SNP groups. Missing columns: \n",paste(missing_cols,collapse = ', '))
    }
   
  }
  
  coloc_dat <- eQTL_Catalogue.annotate_tissues(dat = coloc_dat)  
  # Heatmap
  gg_coloc <- ggplot(data=coloc_dat, aes(x=Locus.eGene, y=qtl_ID, fill=PP.Hyp4)) + 
    # annotate("rect", xmin=-Inf, xmax=Inf, ymin=Tissue_group, ymax=Tissue_group) +
    # geom_rect(data = coloc_dat, aes(xmin=-Inf, xmax=Inf, ymin= qtl_ID, ymax=qtl_ID, color=Tissue_group), alpha=0.4, inherit.aes = F) +
    geom_tile(stat = "identity") +   
    
    # scale_fill_continuous(limits=c(minPP4,1)) +
    scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue",.7), high = scales::alpha("red",.7)) +
    scale_y_discrete(drop=T) +
    scale_x_discrete(position = "top") + 
    facet_grid(facets = System ~ .,
               switch = "y",scales = "free", space = "free") +
    labs(title=paste0("Colocalized GWAS loci x QTL loci (> ",PP_thresh*100,"% probability)\n"),
         x="GWAS Locus (QTL eGene)",
         y="QTL Dataset\n",
         fill="Colocalization\nProbability") +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5),
          axis.text.x = element_text(angle=45, hjust=0), 
          # axis.text.y = element_text(ifelse(unique(coloc_dat$Tissue_group) == 'T-cell', 'red', 'black')),
          panel.border = element_rect(colour = "transparent"),
          strip.placement = "outside", 
          strip.background = element_rect(fill="transparent"), 
          strip.text.y = element_text(angle=90))
  
  if(label_snp_groups){
    gg_coloc <- gg_coloc + 
      # Consensus SNP markers
      geom_point(data = subset(coloc_dat, Consensus.sigQTL>0), aes(x=Locus.eGene, y=qtl_ID, color="Consensus_SNP"), color="cyan2", shape=16, size=1.5, show.legend = F) +
      # UCS SNP markers
      geom_point(data = subset(coloc_dat, UCS.sigQTL>0), aes(x=Locus.eGene, y=qtl_ID, color="Union_Credible_Set"), size=3, shape=5, color="cyan2", show.legend = F) +
      # lead GWAS SNP markers
      geom_tile(data = subset(coloc_dat, leadGWAS.sigQTL>0), aes(x=Locus.eGene, y=qtl_ID), fill="transparent", color="black", size=.7) 
    # geom_point(data = subset(coloc_dat, leadGWAS.sigQTL>0), aes(x=Locus.eGene, y=qtl_ID, color="Lead_GWAS_SNP"), size=5, shape=0, color="black", stroke=1, show.legend=F)
  } 
  if(show_plot) print(gg_coloc) 
  
  
  if(save_dir!=F){
    # save_dir = "./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC"
    save_path <- file.path(save_dir,paste0("coloc_PP",PP_thresh*100,".png"))
    dir.create(dirname(save_path), showWarnings = F, recursive = T)
    if(PP_thresh==.99){
      width <- 7; height <- 9
    } else {
      width <- length(unique(coloc_dat$qtl_ID)) * .23#.333333
      height <- length(unique(coloc_dat$Locus.eGene)) * .3 # 0.6428571
    }
    
    ggsave(save_path,
           plot=gg_coloc, height=width, width=height, dpi=400)
  }
  # Return merged data
  return(gg_coloc)
}




# plot_coloc_summary <- function(coloc_QTLs, 
#                                topQTL,
#                                merged_DT,
#                                PP_thresh = .8, 
#                                save_dir=".",
#                                no_no_loci=NULL,
#                                gwas_dataset="./Data/GWAS/Nalls23andMe_2019",
#                                label_snp_groups=T){  
#   topQTL <- data.table::fread("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_topHits.tsv.gz", nThread = 4)
#   no_no_loci = c("HLA-DRB5","ATG14","SP1","LMNB1","ATP6V0A1", "RETREG3","UBTF","FAM171A2","MAP3K14","CRHR1","MAPT-AS1","KANSL1","NSF","WNT3"); label_snp_groups=T;gwas_dataset="./Data/GWAS/Nalls23andMe_2019";
#   coloc_QTLs <- data.table::fread("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/coloc.eQTL_Catalogue_ALL.csv.gz", nThread = 4)
#   # topQTL <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_sigHits6e-11_finemapOverlap.tsv", nThread = 4)
#   # merged_DT <- echolocatoR::merge_finemapping_results(dataset = gwas_dataset,
#   #                                                     minimum_support = 1,
#   #                                                     include_leadSNPs = T,
#   #                                                     xlsx_path = F)
#   merged_DT <- openxlsx::read.xlsx("~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/Nalls23andMe_2019_results.xlsx")
#   merged_DT <- echolocatoR::update_CS_cols(merged_DT) 
#   
#   # ---------------------------- Add old pipeline results 
#   # coloc_QTLs <- add_old_coloc_results(old_results_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC/old_results/COLOC_results_flip-gwasEffect.txt",
#   #                                     coloc_QTLs)
#   # ---------------------------- 
#   
#   sigQTL <- top_eVariants_overlap(topQTL = topQTL, 
#                                   merged_DT = merged_DT,
#                                   gwas_dataset = gwas_dataset, 
#                                   save_dir = F) 
#   # PP_thresh <- .9
#   
#   coloc_plot <- subset(coloc_QTLs, !is.na(Locus.GWAS)) %>% 
#     dplyr::select(-qtl.ID) %>%
#     dplyr::rename(qtl.ID=qtl.id) %>%
#     dplyr::mutate(PP.H4.thresh = ifelse(PP.H4>=PP_thresh, PP.H4,NA),
#                   PP.Hyp4= ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4/PP.H3 >= 2), PP.H4,NA)) %>%
#     tidyr::separate(col = qtl.ID, into = c("QTL.group","id"), sep = "\\.", remove = F) %>%
#     # only plot colocalized loci-egene combinations (otherwise, plot will be massive)
#     subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop=F) %>% 
#     subset(!Locus.GWAS %in% no_no_loci) %>%
#     data.table::data.table()
#   
#   
#   # if(label_snp_groups){ 
#   coloc_dat <- data.table:::merge.data.table(x=data.table::data.table(coloc_plot),
#                                              y=data.table::data.table(sigQTL),
#                                              by.x=c("Locus.GWAS","eGene","qtl.ID"),
#                                              by.y=c("Locus","eGene","qtl.ID"), 
#                                              all.x=T) %>%
#     dplyr::group_by(Locus.GWAS,eGene,qtl.ID) %>% 
#     dplyr::summarise(PP.Hyp4=max(PP.Hyp4),
#                      leadGWAS.sigQTL = sum(leadSNP, na.rm=T),
#                      Consensus.sigQTL = sum(Consensus_SNP, na.rm=T),
#                      UCS.sigQTL = n_distinct(SNP[Support>0],na.rm = T)
#     ) %>% 
#     subset(!is.na(PP.Hyp4)) %>%
#     data.table::data.table()
#   # coloc_dat
#   # }
#   
#   coloc_dat <- eQTL_Catalogue.annotate_tissues(dat = coloc_dat) 
#   # Do manually for old colocs results 
#   coloc_dat[grep("GTEx_V7.Brain", coloc_dat$qtl.ID),"Tissue"] <- "brain" 
#   coloc_dat[grep("MESA.", coloc_dat$qtl.ID),"Tissue"] <- "monocytes" 
#   coloc_dat[grep("GTEx_V7.Brain", coloc_dat$qtl.ID),"System"] <- "CNS"
#   coloc_dat[grep("MESA.", coloc_dat$qtl.ID),"System"] <- "Blood"
#   
#   coloc_dat$Tissue_count <- gsub(" [(]n=","\n(n=",coloc_dat$Tissue_count)
#   coloc_dat <- coloc_dat %>% dplyr::arrange(System, Tissue)
#   coloc_dat$qtl.ID <- factor(coloc_dat$qtl.ID, unique(coloc_dat$qtl.ID), ordered = T)
#   coloc_dat$Locus.eGene <- paste0(coloc_dat$Locus.GWAS,"  (",coloc_dat$eGene,")")
#   max_x <- length(unique(coloc_dat$Locus.eGene))
#   coloc_dat <- subset(coloc_dat, !is.na(eGene) & eGene!="NA")
#   # colSums(coloc_dat[,c("leadGWAS.sigQTL","Consensus.sigQTL","UCS.sigQTL")])
#   
#   # Heatmap
#   gg_coloc <- ggplot(data=coloc_dat, aes(x=Locus.eGene, y=qtl.ID, fill=PP.Hyp4)) + 
#     # annotate("rect", xmin=-Inf, xmax=Inf, ymin=Tissue_group, ymax=Tissue_group) +
#     # geom_rect(data = coloc_dat, aes(xmin=-Inf, xmax=Inf, ymin= qtl.ID, ymax=qtl.ID, color=Tissue_group), alpha=0.4, inherit.aes = F) +
#     geom_tile(stat = "identity") +   
#     
#     # scale_fill_continuous(limits=c(minPP4,1)) +
#     scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue",.7), high = scales::alpha("red",.7)) +
#     scale_y_discrete(drop=T) +
#     scale_x_discrete(position = "top") + 
#     facet_grid(facets = System ~ .,
#                switch = "y",scales = "free", space = "free") +
#     labs(title=paste0("Colocalized GWAS loci x QTL loci (> ",PP_thresh*100,"% probability)\n"),
#          x="GWAS Locus (QTL eGene)",
#          y="QTL Dataset\n",
#          fill="Colocalization\nProbability") +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = .5),
#           axis.text.x = element_text(angle=45, hjust=0), 
#           # axis.text.y = element_text(ifelse(unique(coloc_dat$Tissue_group) == 'T-cell', 'red', 'black')),
#           panel.border = element_rect(colour = "transparent"),
#           strip.placement = "outside", 
#           strip.background = element_rect(fill="transparent"), 
#           strip.text.y = element_text(angle=90)) +
#     
#     # coord_cartesian(xlim=c(0,max_x+2), clip = 'off') + 
#     # annotate("text", x=max_x+2, y=2.5, label = "Significant eVariant\noverlaps with ≥ 1:", vjust = -1.5, hjust = 0, colour = "black", size=3.5) +
#     # annotate("text", x=max_x+2, y=2.5, label = "•  Consensus SNP", vjust = 2, hjust = 0, colour = "goldenrod2", size=3) +
#     # annotate("text",x=max_x+2, y=2.5, label = "○ Union Credible Set SNP", vjust = 0, hjust = 0, colour = "green3", size=3) +
#     # annotate("text", x=max_x+2, y=2.5, label = "○ Lead GWAS SNP", vjust = -2,  hjust = 0, colour = "black", size=3) + 
#     # Consensus SNP markers
#     geom_point(data = subset(coloc_dat, Consensus.sigQTL>0), aes(x=Locus.eGene, y=qtl.ID, color="Consensus_SNP"), color="cyan2", shape=16, size=1.5, show.legend = F) +
#     # UCS SNP markers
#     geom_point(data = subset(coloc_dat, UCS.sigQTL>0), aes(x=Locus.eGene, y=qtl.ID, color="Union_Credible_Set"), size=3, shape=5, color="cyan2", show.legend = F) +
#     # lead GWAS SNP markers
#     geom_tile(data = subset(coloc_dat, leadGWAS.sigQTL>0), aes(x=Locus.eGene, y=qtl.ID), fill="transparent", color="black", size=.7) 
#   # geom_point(data = subset(coloc_dat, leadGWAS.sigQTL>0), aes(x=Locus.eGene, y=qtl.ID, color="Lead_GWAS_SNP"), size=5, shape=0, color="black", stroke=1, show.legend=F)
#   print(gg_coloc) 
#   
#   
#   if(save_dir!=F){
#     # save_dir = "./Data/GWAS/Nalls23andMe_2019/_genome_wide/COLOC"
#     save_path <- file.path(save_dir,paste0("coloc_PP",PP_thresh*100,".png"))
#     dir.create(dirname(save_path), showWarnings = F, recursive = T)
#     if(PP_thresh==.99){
#       width <- 7; height <- 9
#     } else {
#       width <- length(unique(coloc_dat$qtl.ID)) * .23#.333333
#       height <- length(unique(coloc_dat$Locus.eGene)) * .3 # 0.6428571
#     }
#     
#     ggsave(save_path,
#            plot=gg_coloc, height=width, width=height, dpi=400)
#   }
#   # Return merged data
#   return(coloc_plot)
# }





