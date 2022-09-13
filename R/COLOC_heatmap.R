#' Plot summary of coloc tests: heatmap
#'
#' Plots the output from \link[catalogueR]{COLOC_run}.
#' 
#' @param label_snp_groups Whether to label SNP groups.
#' @param save_dir Directory to save plot in.
#' @param show_plot Whether to print the plot. 
#' @inheritParams COLOC_merge_res
#' 
#' @family coloc
#' @export
#' @importFrom dplyr rename mutate n_distinct group_by summarise
#' @importFrom tidyr separate
#' @importFrom data.table data.table
#' @examples
#' coloc_QTLs <- catalogueR::COLOC_get_example_res()
#' gg_coloc <- catalogueR::COLOC_heatmap(
#'     coloc_QTLs = coloc_QTLs,
#'     coloc_thresh = .5)
COLOC_heatmap <- function(coloc_QTLs,
                          coloc_thresh = .8,
                          label_snp_groups = FALSE, 
                          save_dir = tempdir(),
                          show_plot = TRUE,
                          verbose = TRUE) {
  
  requireNamespace("ggplot2") 
  Locus.eGene <- Consensus.sigQTL <- UCS.sigQTL <- leadGWAS.sigQTL <-
    PP.Hyp4 <- qtl_ID <- NULL;

  #### Prepare data ####
  coloc_dat <- COLOC_plot_data(coloc_QTLs = coloc_QTLs,
                               coloc_thresh = coloc_thresh,
                               label_snp_groups = label_snp_groups, 
                               verbose = verbose)
  #### Heatmap ####
  gg_coloc <- ggplot2::ggplot(
    data = coloc_dat,
    ggplot2::aes(
      x = Locus.eGene,
      y = qtl_ID, fill = PP.Hyp4)
  ) +
    ggplot2::geom_tile(stat = "identity") +
    
    # scale_fill_continuous(limits=c(minPP4,1)) +
    ggplot2::scale_fill_gradient(
      na.value = "transparent",
      low = scales::alpha("blue", .7),
      high = scales::alpha("red", .7)
    ) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::facet_grid(
      facets = System ~ .,
      switch = "y", scales = "free", space = "free"
    ) +
    ggplot2::labs(
      title = paste0(
        "Colocalized GWAS loci x QTL loci (> ",
        coloc_thresh * 100, "% probability)\n"
      ),
      x = "GWAS Locus (QTL eGene)",
      y = "QTL Dataset\n",
      fill = "Colocalization\nProbability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = .5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0),
      panel.border = ggplot2::element_rect(colour = "transparent"),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(fill = "transparent"),
      strip.text.y = ggplot2::element_text(angle = 90)
    )
  
  if (isTRUE(label_snp_groups)) {
    gg_coloc <- gg_coloc +
      # Consensus SNP markers
      ggplot2::geom_point(
        data = subset(coloc_dat, Consensus.sigQTL > 0),
        ggplot2::aes(
          x = Locus.eGene, y = qtl_ID,
          color = "Consensus_SNP"
        ),
        color = "cyan2", shape = 16, size = 1.5,
        show.legend = FALSE
      ) +
      # UCS SNP markers
      ggplot2::geom_point(
        data = subset(coloc_dat, UCS.sigQTL > 0),
        ggplot2::aes(
          x = Locus.eGene, y = qtl_ID,
          color = "Union_Credible_Set"
        ),
        size = 3, shape = 5, color = "cyan2",
        show.legend = FALSE
      ) +
      # lead GWAS SNP markers
      ggplot2::geom_tile(
        data = subset(coloc_dat, leadGWAS.sigQTL > 0),
        ggplot2::aes(x = Locus.eGene, y = qtl_ID),
        fill = "transparent", color = "black", size = .7
      )
  }
  #### Print ####
  if (isTRUE(show_plot)) methods::show(gg_coloc)
  #### Save ####
  save_path <- COLOC_save_plot(gg_coloc=gg_coloc,
                               save_dir=save_dir,
                               coloc_dat=coloc_dat,
                               coloc_thresh=coloc_thresh)
  #### Return merged data ####
  return(list(plot=gg_coloc,
              path=save_path)) 
}

#### Deprecation function #####
plot_coloc_summary <- function(...){
  .Deprecated("COLOC_heatmap")
  COLOC_heatmap(...)
}
