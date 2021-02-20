

#' Get RSIDS from SNP coordinates 
#' 
#' @param dat Data.table of SNP positions, 
#' with the columns \code{CHR} and \code{POS}. Extra columns allowed.
#' @param genome_build Which genome build \code{dat} is in (e.g. "hg19" or "hg38")
#' @param drop_unannotated Drop SNPs that RSIDs couldn't be found for.
#' @param drop_duplicates Drop any duplicate SNPs rows. 
#' @examples 
#' dat <- catalogueR::BST1
#' dat.annot <- rsids_from_coords(dat,  genome_build="hg19")
rsids_from_coords <- function(dat, 
                              genome_build="hg19",
                              drop_unannotated=T,
                              drop_duplicates=T,
                              verbose=T){
  # dat <- catalogueR::BST1
  printer("Searching for RSIDs using",genome_build,".",v=verbose)
  if(tolower(genome_build) %in% c("hg19","grch37") ){
    db <- SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37
  }
  if(tolower(genome_build) %in% c("hg38","grch38")){
    db <- SNPlocs.Hsapiens.dbSNP144.GRCh38::SNPlocs.Hsapiens.dbSNP144.GRCh38
  }
  gr.snp <- GenomicRanges::makeGRangesFromDataFrame(dat ,keep.extra.columns = T, 
                                                    seqnames.field = "CHR", 
                                                    start.field = "POS", 
                                                    end.field = "POS")
  gr.rsids <- BSgenome::snpsByOverlaps(db, ranges = gr.snp, )
  rsids <- data.table::data.table(data.frame(gr.rsids))
  rsids$seqnames <- tolower(as.character(rsids$seqnames))
  dat$CHR <- tolower(as.character(dat$CHR))
  # Merge
  dat.annot <- data.table::merge.data.table(dat, 
                                            rsids, 
                                            all.x = !drop_unannotated,
                                            by.x = c("CHR","POS"),
                                            by.y = c("seqnames","pos")) 
  if(drop_duplicates){
    dat.annot <- dat.annot[!duplicated(dat.annot$RefSNP_id),]
  }
 
  printer(nrow(dat.annot),"/",nrow(dat),"SNPs annotated with RSIDs.",
          v=verbose)
  return(dat.annot)
}



## Old version of function that required installation of motifbreakR
# get_snp_coordinations <- function(snp_list){
#   library(BSgenome)
#   snp_list <- unique(as.character(snp_list[!is.na(snp_list)]))
#   variants <- motifbreakR::snps.from.rsid(rsid = snp_list,
#                                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
#                                           search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19);
#   snp_list
#   return(variants)
# }


