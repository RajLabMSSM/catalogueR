#' Read a tabix file subset
#'
#' Query a subset of a (remote) tabix file and direclty import it with
#' \link[data.table]{fread}.
#'
#' @keywords internal
#' @importFrom echoconda find_package
#' @importFrom data.table fread
tabix_reader <- function(tabix_path,
                         region,
                         conda_env = "echoR",
                         nThread = 1,
                         verbose = TRUE) {
    ## IMPORTANT!! Selecting the appropriate
    ## download.file method now matters for fread.
    ## See details: https://stackoverflow.com/questions/49357652/why-does-it-take-more-time-for-data-tablefread-to-read-a-file-when-filename-is/49357808

    # get_os()=="osx"
    options(download.file.method = "curl")
    tabix <- echoconda::find_package(
        package = "tabix",
        conda_env = conda_env
    )
    tabix_cmd <- paste(
        tabix,
        # "--print-header",
        # "-f",
        tabix_path,
        region
    )
    messager(tabix_cmd, v = verbose)

    # Read directly into R rather than saving tabix subset
    qtl.subset <- data.table::fread(
        cmd = tabix_cmd,
        nThread = nThread,
        header = FALSE
    )
    return(qtl.subset)

    #### Download to tmp then import ####
    # messager("data.table<1.14 detected. Saving tmp file and then reading into R.",v=verbose)
    # token <- generate_token(len=10)
    # tmp_file <- file.path("/tmp/eQTL_Catalogue_queries",paste0(token,".tsv"))
    # dir.create(dirname(tmp_file), recursive  = TRUE, showWarnings  = FALSE)
    # system(paste(tabix_cmd,">",tmp_file))
    # qtl.subset <- data.table::fread(tmp_file,
    #                                 nThread = nThread,
    #                                 header = FALSE)
    #

    #### Rsamtools ####
    ## Rsamtools doesn't know how to parse files properly.
    ## Also tends to freeze R entirely.
    # tbx_con <- Rsamtools::TabixFile(tabix_path)
    # gr.gwas <- GenomicRanges::GRanges(seqnames = chrom,
    #                        ranges =  IRanges::IRanges(bp_lower, bp_upper)  )
    # Rsamtools::scanTabix(tbx_con, gr.gwas)
}
