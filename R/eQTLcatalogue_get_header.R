# Get the header (according to email with Kaur Alassoo)
eQTLcatalogue_get_header <- function(tabix_path=meta$ftp_path[[1]], 
                                     force_new_header=FALSE){
  if(isFALSE(force_new_header)){
    header <- catalogueR::eQTLcatalogue_header
  } else { 
    #### Read in header data #### 
    header <-  colnames(
      data.table::fread(cmd = paste("curl -s",
                                    #### Important! unable to read in using ftp
                                    fix_ftp(tabix_path),
                                    "| zcat | head -n 1"))
      ) 
  }
  return(header)
}