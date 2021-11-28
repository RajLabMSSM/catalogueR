#' Find Consensus SNPs in \emph{echolocatoR} output
#'
#' @family echolocatoR
#' @examples
#' data("BST1")
#' BST1 <- find_consensus_SNPs(finemap_dat = BST1)
find_consensus_SNPs <- function(finemap_dat,
                                verbose = TRUE,
                                credset_thresh = .95,
                                consensus_thresh = 2,
                                sort_by_support = TRUE,
                                exclude_methods = NULL) {
    messager("+ Identifying Consensus SNPs...", v = verbose)
    exclude_methods <- append(exclude_methods, "mean")
    # Find SNPs that are in the credible set for all fine-mapping tools
    CS_cols <- colnames(finemap_dat)[endsWith(colnames(finemap_dat), ".CS")]
    CS_cols <- CS_cols[!(CS_cols %in% paste0(exclude_methods, ".CS"))]
    if (consensus_thresh == "all") {
        consensus_thresh <- length(CS_cols)
    }
    messager("++ support_thresh =", consensus_thresh)
    # Get the number of tools supporting each SNP
    ## Make sure each CS is set to 1
    support_sub <- subset(finemap_dat, select = CS_cols) %>% data.frame()
    support_sub[sapply(support_sub, function(e) {
        e > 1
    })] <- 1
    finemap_dat$Support <- rowSums(support_sub, na.rm = TRUE)
    finemap_dat$Consensus_SNP <- finemap_dat$Support >= consensus_thresh
    # Sort
    if (sort_by_support) {
        finemap_dat <- finemap_dat %>% arrange(desc(Consensus_SNP), desc(Support))
    }

    # Calculate mean PP
    messager("+ Calculating mean Posterior Probability (mean.PP)...")
    PP.cols <- grep(".PP", colnames(finemap_dat), value = TRUE)
    PP.cols <- PP.cols[!(PP.cols %in% paste0(exclude_methods, ".PP"))]
    PP.sub <- subset(finemap_dat, select = c("SNP", PP.cols)) %>% data.frame() # %>% unique()
    PP.sub[is.na(PP.sub)] <- 0
    if (NCOL(PP.sub[, -1]) > 1) {
        finemap_dat$mean.PP <- rowMeans(PP.sub[, -1])
    } else {
        finemap_dat$mean.PP <- PP.sub[, -1]
    }
    finemap_dat$mean.CS <- ifelse(finemap_dat$mean.PP >= credset_thresh, 1, 0)

    # PP.sub %>% arrange(desc(mean.PP)) %>% head()
    messager("++", length(CS_cols), "fine-mapping methods used.")
    messager("++", dim(subset(finemap_dat, Support > 0))[1], "Credible Set SNPs identified.")
    messager("++", dim(subset(finemap_dat, Consensus_SNP = TRUE))[1], "Consensus SNPs identified.")
    return(finemap_dat)
}
