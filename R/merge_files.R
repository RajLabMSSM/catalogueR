#' Merge files from a list of paths
#'
#' Merge a list of files into one by stacking them on top of each other
#'  (i.e. \code{rbind}).
#'  
#' @param file_paths Paths of files to import and merge into one 
#' \link[data.table]{data.table}.
#' @param nThread Number of threads to parallelize reading in files across.
#' @param verbose Print messages. 
#' 
#' @family utils
#' @export
#' @importFrom parallel mclapply
#' @importFrom data.table data.table fread rbindlist
#' @examples
#' sumstats_paths <- echodata:: get_Nalls2019_loci(limit_snps = 5)
#' merged_dat <- merge_files(file_paths = sumstats_paths)
merge_files <- function(file_paths,
                        nThread = 1,
                        verbose = TRUE) {
    messager("+ Merging", length(file_paths), "files.", v = verbose)
    messager("+ Using", nThread,
        paste0(ifelse(nThread > 1, "cores.", "core.")),
        v = verbose
    )
    DAT <- parallel::mclapply(file_paths, function(x) {
        dat <- data.table::data.table()
        try({
            if (endsWith(tolower(x), ".rds")) {
                dat <- readRDS(x)
            } else {
                dat <- data.table::fread(x, nThread = 1)
            }
            dat$file <- basename(x)
        })
        return(dat)
    }, mc.cores = nThread) |> data.table::rbindlist(fill = TRUE)
    messager("+ Merged data.table:", formatC(nrow(DAT), big.mark = ","),
        "rows x", formatC(ncol(DAT), big.mark = ","),
        "columns.",
        v = verbose
    )
    return(DAT)
}

gather_files <- function(...){
    .Deprecated("merge_files")
    merge_files(...)
}
