

# 5. Merge results
make_assocs_list_to_merged_plottable <- function(all_coloc_dt,
                                                 quant_method = "ge"){
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




gather_top_eVariants <- function(root_dir="/pd-omics/brian/eQTL_Catalogue/Nalls23andMe_2019", 
                                 save_path="./eQTL_Catalogue_topHits.tsv.gz",
                                 nThread=4,
                                 criterion="top_eVariant",
                                 pval_thresh=1e-5){ 
  # root_dir = "/Volumes/Steelix/eQTL_Catalogue/Nalls23andMe_2019"; nThread=4;
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
          print(dim(top_qtls))
          
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
# topQTL <- gather_top_eVariants(root_dir="../eQTL_Catalogue/Nalls23andMe_2019", 
#                                             save_path="./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_sigHits.tsv.gz",
#                                             nThread=4,
#                                             criterion="pval_thresh",
#                                             pval_thresh=1e-5)




top_eVariants_overlap <- function(topQTL, 
                                  merged_DT,
                                  gwas_dataset="./Data/GWAS/Nalls23andMe_2019",
                                  gwas_min_support=1,
                                  qtl_pvalue_thresh=NULL, 
                                  save_dir="."){
  # topQTL <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_topHits.tsv.gz", nThread = 4)
  # topQTL <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_sigHits6e-11_finemapOverlap.tsv", nThread = 4)
  
  qtl.cols <- grep(".QTL$",colnames(topQTL), value = T)
  topQTL <- data.table:::merge.data.table(x=data.table::data.table(merged_DT),
                                          y=subset(topQTL, select=c("Locus",'SNP',"qtl.ID", qtl.cols, "eGene")),
                                          by=c("Locus","SNP"),
                                          all=T)
  
  if(is.null(qtl_pvalue_thresh)){
    qtl_pvalue_thresh <- 1e-5 / length(unique(topQTL$qtl.ID)) /length(unique(topQTL$eGene))
  }
  printer("+ Filtering: p <",qtl_pvalue_thresh)
  sigQTL <- topQTL %>%
    subset(!is.na(pvalue.QTL) & pvalue.QTL<qtl_pvalue_thresh) %>%
    # dplyr::group_by(qtl.ID, Locus, eGene) %>% top_n(n=1, wt = -pvalue.QTL) %>%
    dplyr::mutate(Locus.eGene = paste0(Locus,"  (",eGene,")")) %>%
    subset(Support>0 | leadSNP) %>%
    data.table::data.table()
  
  if(save_dir!=F){
    # save_dir="./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue"
    save_path <- file.path(save_dir, paste0("eQTL_Catalogue_sigHits",format(qtl_pvalue_thresh, digits=1),"_finemapOverlap.tsv"))
    data.table::fwrite(qtl_sig, save_path, sep="\t")
  }
  return(sigQTL)
}



top_eVariants_overlap_plot <- function(sigQTL,
                                       locus_order=NULL,
                                       no_no_loci=NULL){ 
  sigQTL <- eQTL_Catalogue.annotate_tissues(sigQTL)
  sigQTL_count <- sigQTL %>% 
    subset(!Locus %in% no_no_loci) %>%  
    dplyr::group_by(Locus, Tissue) %>%  
    dplyr::summarise('Lead GWAS SNPs' = n_distinct(SNP[leadSNP], na.rm=T),
                     "Consensus SNPs" = n_distinct(SNP[Consensus_SNP], na.rm=T),
                     "UCS SNPs" = n_distinct(SNP[Support>0],na.rm = T) ) %>%  
    data.table::data.table() %>%
    data.table::melt.data.table(id.vars = c("Locus","Tissue"), 
                                variable.name = "SNP Group", 
                                value.name = "Significant eVariant Overlap")
  
  sigQTL_count[sigQTL_count$`Significant eVariant Overlap`==0, "Significant eVariant Overlap"] <- NA
  
  if(!is.null(locus_order)){
    sigQTL_count$Locus <- factor(sigQTL_count$Locus,  levels = levels(locus_order$Locus), ordered = T)
  }
  
  ggplot(sigQTL_count, aes(x = Tissue, y = Locus, fill=`Significant eVariant Overlap`)) + 
    geom_tile() + 
    scale_y_discrete(drop=F) + 
    scale_fill_gradient(low = "blue", high="yellow", na.value = "transparent") + 
    facet_grid(facets = .~`SNP Group`) + 
    scale_x_discrete(position = "top") + 
    labs(fill="SNPs overlapping\nwith Significant eVariants") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0)) 
}






file_names <- function(topQTL, root_dir="./"){
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


plot_top_eVariants_overlap <- function(topQTL,
                                       save_path="./topQTL_allLoci.png",
                                       qtl_pvalue_thresh=NULL,
                                       no_no_loci=NULL){ 
  # topQTL <- data.table::fread("./Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_topHits.tsv.gz", nThread = 4)
  # topQTL <- data.table::fread("/pd-omics/brian/Fine_Mapping/Data/GWAS/Nalls23andMe_2019/_genome_wide/eQTL_Catalogue/eQTL_Catalogue_sigHits6e-11_finemapOverlap.tsv", nThread = 4)
  # no_no_loci = c("HLA-DRB5","ATG14","SP1","LMNB1","ATP6V0A1", "CRHR1","MAPT-AS1","KANSL1","NSF","WNT3")
  
  if(is.null(qtl_pvalue_thresh)){
    qtl_pvalue_thresh <- 5e-8 / length(unique(topQTL$qtl.ID)) /length(unique(topQTL$eGene))
  }
  topQTL_sig <- subset(topQTL, (pvalue.QTL<qtl_pvalue_thresh) & Consensus_SNP & (!Locus %in% no_no_loci))
  
  max_x <- length(unique(topQTL_sig$Locus.eGene))
  
  locus.gene_plot <- ggplot(topQTL_sig, aes(x=Locus.eGene, y=qtl.ID)) + 
    # geom_raster(aes(fill=topQTL.prop)) +  
    geom_raster(aes(fill=pvalue.QTL)) + 
    scale_fill_gradient(low = "red", high = "blue", na.value = "transparent") +
    
    # Consensus SNP markers
    # geom_point(data = subset(topQTL_sig, Consensus_SNP), aes(x=Locus.eGene, y=qtl.ID), color="goldenrod2", show.legend = T) +
    # # UCS SNP markers
    # geom_point(data = subset(topQTL_sig, Support>0), aes(x=Locus.eGene, y=qtl.ID), size=3, shape=1, color="green2", show.legend = T) +
    # lead GWAS SNP markers
    geom_point(data = subset(topQTL_sig, leadSNP), aes(x=Locus.eGene, y=qtl.ID), size=3, shape=1, color="black", show.legend=T) +
    coord_cartesian(xlim=c(0,max_x+2), clip = 'off') +
    # annotate("text", x=max_x+2, y=3, label = "• Consensus SNP", vjust = 2, hjust = 0, colour = "goldenrod2", size=3) +
    # annotate("text",x=max_x+2, y=3, label = "○ Union Credible Set", vjust = 0, hjust = 0, colour = "green2", size=3) +
    annotate("text", x=max_x+2, y=3, label = "○ Lead GWAS SNP", vjust = -2,  hjust = 0, colour = "black", size=3) +
    
    scale_y_discrete(drop = F) +
    # facet_grid(facets = .~ Locus.eGene, scales = "free_x") + 
    scale_x_discrete(position = "top") +
    theme_bw() +
    labs(title=paste0("Significant QTL eVariants (p < ",format(qtl_pvalue_thresh, digits=2),") that overlap with Consensus SNPs"),
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





# Insanely large faceted manhattan plots of QTL datasets
manhattan_facets <- function(gwas.qtl_paths,
                             qtl.thresh = 1e-5,
                             gwas_label="GWAS"){
  # gwas.qtl_paths <- list.files("/pd-omics/brian/eQTL_Catalogue/Nalls23andMe_2019", full.names = T, recursive = T)
  # gwas.qtl_paths <- gwas.qtl_paths[1:10]
  PP.H4.thresh=.8
  coloc_QTLs_sig <- coloc_QTLs %>% dplyr::mutate(PP.H4.thresh = ifelse(PP.H4>=PP_thresh, PP.H4,NA),
                                                 PP.Hyp4= ifelse((PP.H3 + PP.H4 >= PP_thresh) & (PP.H4/PP.H3 >= 2), PP.H4,NA)) %>% 
    subset(!is.na(PP.Hyp4))# %>% 
  # subset(Locus.GWAS %in% c("BIN3","MED12L","LRRK2"))
  top_files <- file_names(topQTL = coloc_QTLs_sig,
                          root_dir = "/pd-omics/brian/eQTL_Catalogue/Nalls23andMe_2019") 
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




# #### Check QTL overlap
# 
# ```{r}
# gwas.qtl <- merge_gwas_qtl(gwas_data = merged_DT,
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

