#' Get data via \pkg{piggyback}
#'
#' Download data stored in GitHub Releases using \pkg{piggyback}.
#' 
#' @param fname File name.
#' @param repo GitHub repository \<owner\>/\<name\>.
#' @param storage_dir Data storage directory.
#' @param overwrite Overwrite previously downloaded data.
#' 
#' @keywords internal
get_data <- function(fname,
                     repo = "RajLabMSSM/catalogueR",
                     storage_dir = tempdir(),
                     overwrite = FALSE) {
    requireNamespace("piggyback")
    tmp <- file.path(storage_dir, fname)
    dir.create(storage_dir,showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(tmp)) {
        Sys.setenv("piggyback_cache_duration" = 10)
        piggyback::pb_download(
            file = fname,
            dest = storage_dir,
            repo = repo,
            overwrite = overwrite
        )
    }
    return(tmp)
}
