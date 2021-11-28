# Get the header (according to email with Kaur Alassoo)
tabix_header <- function(tabix_path = file.path(
                             "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv",
                             "Alasoo_2018/ge/Alasoo_2018_ge_macrophage_naive.all.tsv.gz"
                         ),
                         force_new_header = FALSE) {
    if (force_new_header == FALSE) {
        header <- catalogueR::eQTL_Catalogue.header
    } else {
        # Read in header data
        header <- colnames(data.table::fread(
            cmd = paste("curl -s", tabix_path, "| zcat | head -n 1"),
            nThread = 1
        ))
        # Save header data
        # dir.create(dirname(header_path), showWarnings  = FALSE, recursive  = TRUE)
        # data.table::fwrite(list(header), file=header_path, nThread = 1)
    }
    return(header)
}
