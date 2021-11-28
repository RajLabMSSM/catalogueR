prepare_plot_data <- function(coloc_QTLs,
                              locus_name,
                              QTL_gene = NULL,
                              GWAS.QTL = NULL,
                              LD_reference = "1KGphase3",
                              LD_results_dir = "data/COVID19_GWAS",
                              force_new_LD = FALSE) {
    locus_dat <- subset(
        coloc_QTLs,
        Locus.GWAS == locus_name
    ) %>%
        dplyr::rename(Locus = Locus.GWAS) %>%
        dplyr::mutate(
            Mb = POS / 1000000,
            Effect = 1
        )
    if (!is.null(QTL_gene)) locus_dat <- subset(locus_dat, gene.QTL == QTL_gene)

    #### Annotate SNPs by RSIDs  ####
    ## (not provided in original GWAS file)
    library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
    locus.gr <- GenomicRanges::makeGRangesFromDataFrame(locus_dat,
        keep.extra.columns = TRUE,
        seqnames.field = "CHR",
        start.field = "POS", end.field = "POS"
    )
    rsids <- BSgenome::snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh38, locus.gr)
    locus_dat <- merge(locus_dat,
        data.frame(rsids) %>%
            dplyr::rename(SNP = RefSNP_id) %>%
            dplyr::mutate(seqnames = as.integer(seqnames)),
        by.x = c("CHR", "POS"), by.y = c("seqnames", "pos"), all.x = TRUE
    )

    #### Liftover ####
    locus_lift <- catalogueR::liftover(
        gwas_data = locus_dat,
        build.conversion = "hg38.to.hg19"
    ) %>%
        data.frame() %>%
        dplyr::rename(POS = start)
    if (!"leadSNP" %in% colnames(locus_lift)) locus_lift <- echolocatoR::assign_lead_SNP(locus_lift)


    LD_out <- echolocatoR::LD.load_or_create(
        locus_dir = file.path(LD_results_dir, locus_lift$Locus[1]),
        force_new_LD = force_new_LD,
        subset_DT = locus_lift,
        LD_reference = LD_reference
    )
    LD_matrix <- LD_out$LD
    # locus_lift <- LD_out$DT
    plot_dat <- echolocatoR::LD.get_lead_r2(
        finemap_dat = locus_lift,
        LD_matrix = as.matrix(LD_matrix),
        LD_format = "matrix"
    ) %>%
        catalogueR::eQTL_Catalogue.annotate_tissues() %>%
        tidyr::separate(
            col = "qtl_id", sep = "[.]", remove = FALSE,
            into = c("qtl", "tissue_condition"), extra = "drop"
        )
    if (!is.null(GWAS.QTL)) {
        # Merge QTL FDR into coloc results
        sig_qtls <- GWAS.QTL %>%
            subset(qtl_id %in% unique(plot_dat$qtl_id) &
                gene.QTL %in% unique(plot_dat$gene.QTL) &
                FDR.QTL < .05) %>%
            tidyr::separate(
                col = "qtl_id", sep = "[.]", remove = FALSE,
                into = c("qtl", "tissue_condition"), extra = "drop"
            ) %>%
            dplyr::mutate(
                Mb = POS / 1000000,
                sig_QTL = TRUE
            )
        plot_dat <- merge(plot_dat,
            sig_qtls[, c("SNP", "qtl_id", "gene.QTL", "FDR.QTL", "sig_QTL")],
            all.x = TRUE,
            by = c("SNP", "qtl_id", "gene.QTL")
        )
    }
    return(plot_dat)
}
