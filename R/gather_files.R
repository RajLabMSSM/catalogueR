#' Merge files from a list of paths
#'
#' Merge a list of files into one by stacking them on top of each other (i.e. rbind).
#' @family utils
#' @examples
#' sumstats_paths <- example_sumstats_paths()
#' merged_dat <- gather_files(file_paths = sumstats_paths)
gather_files <- function(file_paths,
                         nThread = 1,
                         verbose = TRUE) {
    messager("+ Merging", length(file_paths), "files.", v = verbose)
    messager("+ Using", nThread, 
             paste0(ifelse(nThread > 1, "cores.", "core.")),
             v = verbose)
    DAT <- parallel::mclapply(file_paths, function(x) {
        dat <- data.table::data.table()
        try({
            if (endsWith(tolower(x), ".rds")) {
                dat <- readRDS(x)
            } else {
                dat <- data.table::fread(x)
            }
            dat$file <- basename(x)
        })
        return(dat)
    }, mc.cores = nThread) %>% data.table::rbindlist(fill = TRUE)
    messager("+ Merged data.table:", nrow(DAT), "rows x", ncol(DAT), "columns.", v = verbose)
    return(DAT)
}
