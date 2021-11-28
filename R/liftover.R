#' Lift genome across builds
#'
#' @param build_conversion "hg19.to.hg38" (\emph{default}) or "hg38.to.hg19.
#' @family utils
#' @examples
#' data("BST1")
#' gr.lifted <- liftover(gwas_data = BST1, build.conversion = "hg19.to.hg38")
liftover <- function(gwas_data,
                     build.conversion = "hg19.to.hg38",
                     verbose = TRUE) {
    messager("XGR:: Lifting genome build:", build.conversion, v = verbose)
    # Save original coordinates and SNP IDs
    gwas_data <- gwas_data %>% dplyr::mutate(
        chrom = paste0("chr", gsub("chr", "", CHR)),
        POS.orig = POS,
        SNP.orig = SNP
    )
    # chain <- rtracklayer::import.chain(con = chain_paths$hg19_to_hg38)
    gr.gwas <- GenomicRanges::makeGRangesFromDataFrame(
        df = gwas_data,
        keep.extra.columns = TRUE,
        seqnames.field = "chrom",
        start.field = "POS",
        end.field = "POS"
    )
    gr.lifted <- xLiftOver(
        data.file = gr.gwas, # dplyr::select(gwas_data, CHR, POS),
        format.file = "GRanges",
        build.conversion = build.conversion,
        verbose = verbose,
        merged = FALSE
    ) # merge must =F in order to work
    return(gr.lifted)
}
