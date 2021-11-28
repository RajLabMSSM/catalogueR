# Get the size of the largest file in a dir
max_file_size <- function(top_dir, recursive = TRUE) {
    #
    files <- list.files(top_dir, full.names = TRUE, recursive = recursive)
    inf <- file.info(files)
    print(paste("File size info (Mb) for", length(files), "files."))
    summary(inf$size / 1e6)
}
