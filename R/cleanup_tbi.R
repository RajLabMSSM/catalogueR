# Remove extra .tbi files
cleanup_tbi <- function(DIR = "./",
                        recursive = FALSE) {
    tbi_files <- list.files(DIR, pattern = ".gz.tbi$", recursive = recursive)
    if (length(tbi_files) > 0) out <- file.remove(tbi_files)
}
