check_tabix_method <- function(method){
  method <- tolower(method)[1] 
  opts <- c("echotabix","tabix")
  if(!method %in% opts){
    stop_msg <- paste0("method must be one of:\n",
                      paste(" -",opts, collapse = "\n"))
    stop(stop_msg)
  } else {
    return(method)
  }
}