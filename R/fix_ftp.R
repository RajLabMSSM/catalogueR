fix_ftp <- function(ftp_path){
  gsub("^ftp:","http:",ftp_path)
}