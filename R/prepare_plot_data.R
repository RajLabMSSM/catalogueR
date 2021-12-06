#' Prepare plot data
#' 
#' Prepare plot data.
#' @keywords internal
#' @importFrom echotabix liftover 
prepare_plot_data <- function(coloc_QTLs,
                              locus_name,
                              QTL_gene = NULL,
                              GWAS.QTL = NULL,
                              LD_reference = "1KGphase3",
                              LD_results_dir = file.path(tempdir(),
                                                         "plot_data"),
                              force_new_LD = FALSE,
                              verbose = TRUE) {
    requireNamespace("echoLD")
    requireNamespace("echodata")
    requireNamespace("dplyr")  
    
    Locus.GWAS <- POS <- BP <- gene.QTL <- RefSNP_id <- seqnames <- start <-
        qtl_id <- FDR.QTL <- NULL;
    locus_dat <- subset(
        coloc_QTLs,
        Locus.GWAS == locus_name
    ) %>%
        dplyr::rename(Locus = Locus.GWAS) %>%
        dplyr::mutate(
            Mb = POS / 1000000,
            Effect = 1
        )
    if (!is.null(QTL_gene)) {
        locus_dat <- subset(locus_dat, gene.QTL == QTL_gene)
    }
    #### Annotate SNPs with RSIDs  ####
    ## (not provided in original GWAS file)
    if(!"SNP" %in% colnames(locus_dat)){
        locus_dat <- echodata::coords_to_rsids(dat = locus_dat,
                                               genome_build = "HG38",
                                               snp_colname = "SNP",
                                               verbose = verbose)
    }
    #### Liftover ####
    locus_lift <- echotabix::liftover( 
        sumstats_dt = locus_dat,
        ref_genome = "HG38",
        convert_ref_genome = "HG19",
        verbose = verbose
    ) %>%  dplyr::rename(POS = BP)
    if (!"leadSNP" %in% colnames(locus_lift)) {
        locus_lift <- echodata::assign_lead_SNP(locus_lift)
    }

    LD_out <- echoLD::load_or_create(
        locus_dir = file.path(LD_results_dir, locus_lift$Locus[1]),
        force_new_LD = force_new_LD,
        dat = locus_lift,
        LD_reference = LD_reference
    )
    LD_matrix <- LD_out$LD
    plot_dat <- echoLD:::get_lead_r2(
        dat = locus_lift,
        LD_matrix = as.matrix(LD_matrix),
        LD_format = "matrix"
    ) %>%
        eQTL_Catalogue.annotate_tissues() %>%
        tidyr::separate(
            col = "qtl_id", sep = "[.]", remove = FALSE,
            into = c("qtl", "tissue_condition"), extra = "drop"
        )
    if (!is.null(GWAS.QTL)) {
        #### Merge QTL FDR into coloc results ####
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
