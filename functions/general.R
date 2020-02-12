## ---------------- General Functions ----------------  ##

library(EnsDb.Hsapiens.v75)

printer <- function(..., v=T){if(v){print(paste(...))}}

hgnc_to_ensembl <- function(gene_symbols){
  # columns(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = gene_symbols,
                                      keytype = "SYMBOL",
                                      column = "GENEID")
  return(conversion)
}
ensembl_to_hgnc <- function(ensembl_ids){
  ensembl_ids[is.na(ensembl_ids)] <- "NA"
  conversion <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                                      keys = ensembl_ids,
                                      keytype = "GENEID",
                                      column = "SYMBOL")
  return(conversion)
}

createDT <- function(DF, caption="", scrollY=400){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}

createDT_html <- function(DF, caption="", scrollY=400){
  htmltools::tagList( createDT(DF, caption, scrollY))
}



