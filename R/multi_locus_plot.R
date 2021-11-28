#' Faceted Manhattan plots of QTL datasets
#'
#' Use the coloc results (coloc_QTLs) to determine which full summary stats (gwas.qtl_paths) to plot.
#' @family plots
#' @examples
#' data("coloc_QTLs")
#' gwas.qtl_paths <- example_eQTL_Catalogue_query_paths()
#' gg_gwas.qtl <- multi_locus_plot(gwas.qtl_paths = gwas.qtl_paths, coloc_QTLs = coloc_QTLs, coloc_thresh = .5, qtl_thresh = .005, remove_extra_panes = FALSE)
#' \dontrun{
#' data("coloc_QTLs_full")
#' library(ggplot2)
#' library(dplyr)
#' root_dir <- "/pd-omics/brian/catalogueR_queries/Nalls23andMe_2019"
#' gwas.qtl_paths <- list.files(root_dir, pattern = "*.tsv.gz", recursive = TRUE, full.names = TRUE)
#' gwas.qtl_paths <- gsub("___", "__", gwas.qtl_paths)
#' colnames(coloc_QTLs_full) <- gsub("\\.df1", ".QTL", colnames(coloc_QTLs_full))
#' colnames(coloc_QTLs_full) <- gsub("\\.df2", ".GWAS", colnames(coloc_QTLs_full))
#' plot_dat <- gather_colocalized_data(gwas.qtl_paths = gwas.qtl_paths, coloc_QTLs = coloc_QTLs_full, qtl_thresh = NULL, coloc_thresh = .95)
#' gg_gwas.qtl <- multi_locus_plot(plot_dat = plot_dat)
#' }
multi_locus_plot <- function(gwas.qtl_paths = NULL,
                             coloc_QTLs = NULL,
                             plot_dat = NULL,
                             qtl_thresh = 1e-5,
                             coloc_thresh = .8,
                             gwas_label = "GWAS",
                             remove_extra_panes = TRUE,
                             y_facet_angle = 0, # 270,
                             x_facet_angle = 270,
                             show_plot = TRUE,
                             verbose = TRUE) {
    if (is.null(plot_dat)) {
        plot_dat <- gather_colocalized_data(
            gwas.qtl_paths = gwas.qtl_paths,
            coloc_QTLs = coloc_QTLs,
            qtl_thresh = qtl_thresh,
            coloc_thresh = coloc_thresh,
            gwas_label = gwas_label,
            remove_extra_panes = remove_extra_panes,
            verbose = verbose
        )
    }
    # Plot
    alpha <- .2
    gg.gwas <- ggplot(plot_dat, aes(x = MB, y = -log10(P), color = Locus.GWAS)) +
        geom_hline(yintercept = -log10(5e-8), alpha = .5, linetype = "dashed", size = .5) +
        geom_point(size = .5, alpha = alpha) +
        labs(x = NULL) +
        facet_grid(
            facets = GWAS.label ~ Locus.GWAS + gene.QTL, # gsub(paste0("__",gene.QTL)[1],"",Locus__eGene) ,
            scales = "free"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            strip.text.y = element_text(angle = y_facet_angle),
            strip.text.x = element_text(angle = x_facet_angle),
            panel.grid.minor = element_blank()
        )

    ann_text <- subset(plot_dat, !is.na(PP.H4))[, c("Locus.GWAS", "qtl_id", "gene.QTL", "PP.H4")] %>% unique()
    gg.qtls <- ggplot(plot_dat, aes(x = MB, y = -log10(pvalue.QTL))) +
        geom_point(aes(color = PP.H4), alpha = alpha, size = .5, show.legend = TRUE) +
        scale_color_gradient(low = "blue", high = "red", na.value = "grey") +
        geom_hline(yintercept = -log10(qtl_thresh), alpha = .5, linetype = "dashed", size = .5) +
        facet_grid(
            facets = gsub("\\.", "\n", qtl_id) ~ Locus.GWAS + gene.QTL, # gsub(paste0(Locus.GWAS,"__")[1],"",Locus__eGene),
            scales = "free"
        ) +
        geom_text(data = ann_text, aes(x = Inf, y = Inf, label = paste0(round(PP.H4, 2))), hjust = 1.2, vjust = 1.5, alpha = .7) +
        theme_bw() +
        theme(
            strip.background.x = element_rect(fill = "transparent", colour = "transparent"),
            strip.text.x = element_text(angle = x_facet_angle, size = 0, color = "transparent"),
            strip.text.y = element_text(angle = y_facet_angle),
            axis.text.x = element_blank(),
            panel.grid.minor = element_blank()
        )
    # gg.qtls

    # Merge plots
    gg_gwas.qtl <- gg.gwas + gg.qtls +
        patchwork::plot_layout(ncol = 1) # , heights = c(.1,1))
    if (show_plot) print(gg_gwas.qtl)

    return(gg_gwas.qtl)
}










# merged_plot <- function(GWAS.QTL){
#   # Group and melt
#   GWAS.QTL <- find_consensus_SNPs(GWAS.QTL)
#   gwas.qtl.melt <-  data.table::setDT(GWAS.QTL)[, .(Count = data.table::uniqueN(SNP[Support<0],na.rm  = TRUE)),
#                                                 by=c("Locus","qtl_id","SNP","gene.QTL")] #"leadSNP","Consensus_SNP","UCS"
#   gwas.qtl.melt <- gwas.qtl.melt %>%
#     dplyr::group_by(Locus, qtl_id, SNP) %>%
#     dplyr::summarise(Count = sum(Count, na.rm = TRUE),
#                      Genes = paste(gene.QTL,collapse=", "))
#
#   gwas.qtl.melt[gwas.qtl.melt$Count==0,"Count"] <- NA
#   gg_gwas.qtl <- ggplot(data=gwas.qtl.melt, aes(x=qtl_id, y=Locus, fill=Count)) +
#     geom_raster() +
#     # scale_fill_manual(values = consensus_colors) +
#     scale_fill_gradient(na.value = "transparent", low = scales::alpha("blue",.7), high = scales::alpha("red",.7)) +
#     geom_point(aes(size=ifelse(Count>0, "dot", "no_dot")), show.legend  = FALSE, alpha=.8, color="white") +
#     scale_size_manual(values=c(dot=.5, no_dot=NA), guide="none") +
#     labs(fill = "UCS SNP Count") +
#     theme_classic() +
#     theme(legend.position = "top",
#           legend.title.align = .5,
#           axis.text.x = element_text(angle = 40, hjust = 1),
#           # legend.background =  element_rect(fill = "lightgray"),
#           legend.key = element_rect(colour = "gray60"),
#           legend.text = element_text(size = 8),
#           legend.text.align = .5,
#           # legend.key.size = unit(.5, units = "cm" ),
#           legend.box="horizontal",
#           # panel.background = element_rect(fill = 'transparent'),
#           # panel.grid = element_line(color="gray", size=5),
#           panel.grid.major = element_line(color="grey", size=5) ) +
#     guides(color = guide_legend(nrow = 1, reverse  = FALSE,
#                                 title.position = "top",
#                                 # label.position = "top",
#                                 title.hjust = .5,
#                                 label.hjust = -1)) +
#     # Keep unused levels/Loci
#     scale_y_discrete(drop=FALSE)
#   print(gg_gwas.qtl)
# }
