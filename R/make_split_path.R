# Consistent means of making split path
make_split_path <- function(output_dir,
                            dataset_id,
                            loc) {
    split_path <- file.path(
        output_dir,
        dataset_id,
        paste0(loc, "__", dataset_id, ".tsv.gz")
    )
    names(split_path) <- paste(loc, dataset_id, sep = "__")
    dir.create(dirname(split_path), showWarnings = FALSE, recursive = TRUE)
    return(split_path)
}
