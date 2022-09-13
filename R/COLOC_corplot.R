#' Plot summary of coloc tests: heatmap
#'
#' Plots the output from \link[catalogueR]{COLOC_run}.
#' 
#' @param label_top_snps Label the top n SNPs per \code{facets} group, sorted
#' by \code{SNP.PP.H4} from \pkg{coloc}. 
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
#' gg_coloc <- COLOC_corplot(coloc_QTLs = coloc_QTLs,
#'                           coloc_thresh = .5)
COLOC_corplot <- function(coloc_QTLs,
                          coloc_thresh = .8,
                          label_top_snps = 1,
                          facets = Study + System + Tissue  ~  Locus.eGene,
                          show_plot = TRUE,
                          save_dir = tempdir(),
                          seed = 2022,
                          verbose=TRUE) {
  requireNamespace("ggplot2") 
  SNP.PP.H4 <- NULL;
  
  #### Prepare data ####
  coloc_dat <- COLOC_plot_data(coloc_QTLs = coloc_QTLs,
                               coloc_thresh = 0,
                               label_snp_groups = FALSE,
                               verbose = verbose)  
  #### Plot #### 
  gg_coloc <-  ggplot2::ggplot(coloc_dat, 
                  ggplot2::aes(x=-log10(P), 
                               y=-log10(pvalue.QTL),
                               label=SNP)) + 
    ggplot2::geom_density_2d_filled(contour_var = "ndensity") +
    ggplot2::geom_point(ggplot2::aes(size=SNP.PP.H4, 
                                     color=SNP.PP.H4),
                        # shape=16,
                        alpha=0.1) +
    ggplot2::scale_color_gradient(low="darkblue", high = "cyan2",
                                   limits=c(0,1), breaks=c(0,0.5,1)) + 
    ggplot2::scale_fill_viridis_d(option = "magma",) +
    ggplot2::facet_grid(facets = facets,
                        scales = "free") + 
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(color = "black"))
  
   
  if (label_top_snps>0) {
    #### Get top SNPs ####
    top_snps <- 
      coloc_dat |>
      subset(PP.Hyp4>coloc_thresh & SNP.PP.H4>coloc_thresh) |>
      dplyr::group_by_at(.vars = all.vars(facets)) |>
      dplyr::slice_max(order_by = SNP.PP.H4, 
                       n = 1, 
                       with_ties = FALSE) |>
      data.table::data.table()
    if(nrow(top_snps)>0){
      gg_coloc <- gg_coloc + 
        ggrepel::geom_label_repel(data = top_snps,
                                  segment.color="white",
                                  seed = seed,
                                  alpha=.75) 
    } 
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
