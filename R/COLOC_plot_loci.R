#' Faceted Manhattan plots of QTL datasets
#'
#' Use the coloc results (\code{coloc_QTLs}) to determine which
#' full summary statistics (\code{gwas.qtl_paths}) to plot.
#' 
#' @param plot_dat [Optional] Pre-computed plot data from
#'  \link[catalogueR]{COLOC_merge_res}.
#' @param y_facet_angle Angle of the y-axis facet labels.
#' @param x_facet_angle Angle of the x-axis facet labels.
#' @param show_plot Whether to print the plot.
#' @inheritParams COLOC_merge_res
#' 
#' @family plots
#' @export
#' @importFrom patchwork plot_layout
#' @examples
#' coloc_QTLs <- catalogueR::COLOC_get_example_res()
#' gwas.qtl_paths <- catalogueR::eQTLcatalogue_example_queries()
#' gg_gwas.qtl <- catalogueR::COLOC_plot_loci(
#'     gwas.qtl_paths = gwas.qtl_paths,
#'     coloc_QTLs = coloc_QTLs,
#'     coloc_thresh = .5,
#'     qtl_thresh = .005,
#'     remove_extra_panes = FALSE)  
COLOC_plot_loci <- function(gwas.qtl_paths = NULL,
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
    requireNamespace("ggplot2")
    
    MB <- P <- Locus.GWAS <- PP.H4 <- pvalue.QTL <- NULL;
    if (is.null(plot_dat)) {
        plot_dat <- COLOC_merge_res(
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
    gg.gwas <- ggplot2::ggplot(plot_dat, 
                               ggplot2::aes(x = MB, y = -log10(P),
                                            color = Locus.GWAS)) +
        ggplot2::geom_hline(yintercept = -log10(5e-8), alpha = .5, 
                            linetype = "dashed", size = .5) +
        ggplot2::geom_point(size = .5, alpha = alpha) +
        ggplot2::labs(x = NULL) +
        ggplot2::facet_grid(
            facets = GWAS.label ~ Locus.GWAS + gene.QTL,  
            scales = "free"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            strip.text.y = ggplot2::element_text(angle = y_facet_angle),
            strip.text.x = ggplot2::element_text(angle = x_facet_angle),
            panel.grid.minor = ggplot2::element_blank()
        )

    ann_text <- subset(plot_dat, 
                       !is.na(PP.H4))[, c("Locus.GWAS", "qtl_id", 
                                          "gene.QTL", "PP.H4")] |> unique()
    gg.qtls <- ggplot2::ggplot(plot_dat, 
                               ggplot2::aes(x = MB, y = -log10(pvalue.QTL))) +
        ggplot2::geom_point(ggplot2::aes(color = PP.H4),
                            alpha = alpha, size = .5,
                   show.legend = TRUE) +
        ggplot2::scale_color_gradient(low = "blue", high = "red"
                                      , na.value = "grey") +
        ggplot2::geom_hline(yintercept = -log10(qtl_thresh), alpha = .5,
                   linetype = "dashed", size = .5) +
        ggplot2::facet_grid(
            facets = gsub("\\.", "\n", qtl_id) ~ Locus.GWAS + gene.QTL, 
            scales = "free"
        ) +
        ggplot2::geom_text(data = ann_text, 
                  ggplot2::aes(x = Inf, y = Inf, 
                               label = paste0(round(PP.H4, 2))),
                  hjust = 1.2, vjust = 1.5, alpha = .7) +
        ggplot2:: theme_bw() +
        ggplot2:: theme(
            strip.background.x = ggplot2::element_rect(fill = "transparent", 
                                                       colour = "transparent"),
            strip.text.x = ggplot2::element_text(angle = x_facet_angle, 
                                                 size = 0, 
                                                 color = "transparent"),
            strip.text.y = ggplot2::element_text(angle = y_facet_angle),
            axis.text.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )
    # gg.qtls

    # Merge plots
    gg_gwas.qtl <- gg.gwas + gg.qtls +
        patchwork::plot_layout(ncol = 1) # , heights = c(.1,1))
    if (show_plot) print(gg_gwas.qtl)

    return(gg_gwas.qtl)
}


#### Deprecation function #####
multi_locus_plot <- function(...){
  .Deprecated("COLOC_plot_loci")
  COLOC_plot_loci(...)
}


# merged_plot <- function(GWAS.QTL){
#   # Group and melt
#   GWAS.QTL <- find_consensus_SNPs(GWAS.QTL)
#   gwas.qtl.melt <-  data.table::setDT(GWAS.QTL)[, .(Count = data.table::uniqueN(SNP[Support<0],na.rm  = TRUE)),
#                                                 by=c("Locus","qtl_id","SNP","gene.QTL")] #"leadSNP","Consensus_SNP","UCS"
#   gwas.qtl.melt <- gwas.qtl.melt |>
#     dplyr::group_by(Locus, qtl_id, SNP) |>
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
