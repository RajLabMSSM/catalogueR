COLOC_save_plot <- function(gg_coloc,
                            save_dir,
                            coloc_dat,
                            coloc_thresh){
  if (!is.null(save_dir)) {
    save_path <- file.path(save_dir, paste0(
      "coloc_PP",
      coloc_thresh, ".png"
    ))
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    if (coloc_thresh == .99) {
      width <- 7
      height <- 9
    } else {
      width <- length(unique(coloc_dat$qtl_ID)) * .23 # .333333
      height <- length(unique(coloc_dat$Locus.eGene)) * .3 # 0.6428571
    } 
    ggplot2::ggsave(
      filename = save_path,
      plot = gg_coloc,
      height = width, 
      width = height
    )
    return(save_path)
  } else {
    return(NULL)
  }
}
