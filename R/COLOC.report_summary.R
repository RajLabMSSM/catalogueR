
#' Report coloc results
#'
#' @family coloc
#' @keywords internal
COLOC.report_summary <- function(coloc.res,
                                 coloc_thresh = .8) {
    # MAF = dataset1$MAF)
    hypothesis_key <- setNames(
        c(
            "Neither trait has a genetic association in the region.",
            "Only trait 1 has a genetic association in the region.",
            "Only trait 2 has a genetic association in the region.",
            "Both traits are associated, but with different causal variants",
            "(one in each dataset).",
            "Both traits are associated and share a single causal variant."
        ),
        c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    )
    # Report hypothess results
    messager("Hypothesis Results @ coloc_thresh =", coloc_thresh, ":")
    true_hyp <- ""
    for (h in names(hypothesis_key)) {
        if (!is.na(coloc.res$summary[h])) {
            if (coloc.res$summary[h] >= coloc_thresh) {
                hyp <- hypothesis_key[h]
                messager("    ", h, "== TRUE: **", hyp)
                true_hyp <- paste0(names(hyp), ": ", hyp)
            }
        } else {
            messager("    ", h, "== FALSE: ")
        }
    }

    # Save raw results
    coloc_DT <- coloc.res$results
    # Process results
    coloc_DT$Colocalized <- coloc_DT$SNP.PP.H4 >= coloc_thresh
    colocalized_snps <- subset(coloc_DT, Colocalized = TRUE)$snp  
    subtitle2 <- paste0("Colocalized SNPs: ",
                        paste(colocalized_snps, sep = ", "))
    if (!is.na(coloc.res$summary)["PP.H4.abf"]) {
        if ((coloc.res$summary["PP.H3.abf"] + 
             coloc.res$summary["PP.H4.abf"] >= coloc_thresh) &
            (coloc.res$summary["PP.H4.abf"] / 
             coloc.res$summary["PP.H3.abf"] >= 2)) { 
            report <- paste("Datasets colocalized")
        } else {
            report <- paste("Datasets NOT colocalized")
        }
    } else {
        report <- paste("Datasets NOT colocalized")
    }
    messager("+ COLOC::", report, "at: PP.H3 + PP.H4 >=", 
             coloc_thresh, " and PP.H3 / PP.H4 >= 2.")
    return(coloc_DT)
}
