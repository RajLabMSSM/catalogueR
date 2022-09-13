COLOC_plot_data <- function(coloc_QTLs,
                            coloc_thresh,
                            label_snp_groups,
                            verbose=TRUE){
  
  Locus.GWAS <- qtl_id <- PP.H4 <- PP.H3 <- gene.QTL <- qtl_ID <- PP.Hyp4 <- 
    PP.H4.thresh <- eGene <- leadSNP <- Consensus_SNP <- SNP <- Support <- NULL;
  
  messager("Perparing COLOC plot data.",v=verbose)
  if(("snp" %in% names(coloc_QTLs)) && (!"SNP" %in% names(coloc_QTLs))){
    data.table::setnames(coloc_QTLs, "snp","SNP")
  }
  coloc_dat <- subset(coloc_QTLs, !is.na(Locus.GWAS)) |>
    dplyr::rename(qtl_ID = qtl_id) |>
    dplyr::mutate(
      PP.H4.thresh = ifelse(PP.H4 >= coloc_thresh, PP.H4, NA),
      PP.Hyp4 = ifelse((PP.H3 + PP.H4 >= coloc_thresh) &
                         (PP.H4 / PP.H3 >= 2), PP.H4, NA),
      Locus.eGene = paste0(Locus.GWAS, " (", gene.QTL, ")")
    ) |>
    # tidyr::separate(
    #   col = qtl_ID, 
    #   into = c("group.QTL", "id"),
    #   sep = "\\.", remove = FALSE) |>
    # only plot colocalized loci-egene combinations
    # (otherwise, plot will be massive)
    subset((!is.na(PP.Hyp4) & !is.na(PP.H4.thresh)), .drop = FALSE) |>
    data.table::data.table()
  if (nrow(coloc_dat) == 0) {
    stop_msg <- paste(
      "\n\n** No coloc tests reached the PPH4 threshold",
      "(", coloc_thresh, ").",
      "Try lowering `coloc_thresh`. **"
    )
    stop(stop_msg)
  }
  if(!"leadSNP" %in% names(coloc_dat)){
    coloc_dat <- echodata::assign_lead_snp(dat = coloc_dat, 
                                           grouping_vars = "Locus.GWAS",
                                           verbose = verbose)
  }  
  #### Check needed cols ####
  needed_cols <- c("leadSNP", "Consensus_SNP", "Support")
  if (isTRUE(label_snp_groups)) {
    if (all(needed_cols %in% colnames(coloc_dat))) {
      coloc_dat <- coloc_dat |>
        dplyr::group_by(Locus.GWAS, gene.QTL, qtl_ID) |>
        dplyr::summarise(
          PP.Hyp4 = max(PP.Hyp4),
          leadGWAS.sigQTL = sum(leadSNP, na.rm = TRUE),
          Consensus.sigQTL = sum(Consensus_SNP, na.rm = TRUE),
          UCS.sigQTL = dplyr::n_distinct(SNP[Support > 0],
                                         na.rm = TRUE
          )
        ) |>
        subset(!is.na(PP.Hyp4)) |>
        data.table::data.table()
    } else {
      missing_cols <- needed_cols[!needed_cols %in% colnames(coloc_dat)]
      messager(
        "++ Cannot label SNP groups. Missing columns: \n",
        paste(missing_cols, collapse = ", "),v=verbose
      )
    }
  } 
  coloc_dat <- eQTLcatalogue_annotate_tissues(dat = coloc_dat)
  data.table::setkey(coloc_dat,PP.Hyp4)
  return(coloc_dat)
}
