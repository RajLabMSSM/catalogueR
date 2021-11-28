parse_gwas.qtl_path <- function(gwas.qtl_path,
                                get_locus = FALSE,
                                get_qtl_id = TRUE,
                                sep = "__") {
    split_name <- strsplit(basename(gwas.qtl_path), sep)[[1]]
    if (get_locus & get_qtl_id == FALSE) {
        return(split_name[1])
    }
    if (get_locus == FALSE & get_qtl_id) {
        return(gsub(".tsv|.gz", "", split_name[2]))
    }
    if (get_locus & get_qtl_id) {
        return(list(
            locus = split_name[1],
            qtl_id = gsub(".tsv|.gz", "", split_name[2])
        ))
    }
}
