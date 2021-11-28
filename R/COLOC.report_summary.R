
#' Report coloc results
#'
#' @family coloc
#' @keywords internal
COLOC.report_summary <- function(coloc.res,
                                 PP_threshold = .8) {
    # MAF = dataset1$MAF)
    hypothesis_key <- setNames(
        c(
            "Neither trait has a genetic association in the region.",
            "Only trait 1 has a genetic association in the region.",
            "Only trait 2 has a genetic association in the region.",
            "Both traits are associated, but with different causal variants (one in each dataset).",
            "Both traits are associated and share a single causal variant."
        ),
        c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    )
    # Report hypothess results
    messager("Hypothesis Results @ PP_threshold =", PP_threshold, ":")
    true_hyp <- ""
    for (h in names(hypothesis_key)) {
        if (!is.na(coloc.res$summary[h])) {
            if (coloc.res$summary[h] >= PP_threshold) {
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
    coloc_DT$Colocalized <- ifelse(coloc_DT$SNP.PP.H4 >= PP_threshold, T, F)
    colocalized_snps <- subset(coloc_DT, Colocalized = TRUE)$snp # subset(coloc_DT, Colocalized==1)$SNP
    subtitle2 <- paste0("Colocalized SNPs: ", paste(colocalized_snps, sep = ", "))
    if (!is.na(coloc.res$summary)["PP.H4.abf"]) {
        if ((coloc.res$summary["PP.H3.abf"] + coloc.res$summary["PP.H4.abf"] >= PP_threshold) &
            (coloc.res$summary["PP.H4.abf"] / coloc.res$summary["PP.H3.abf"] >= 2)) {
            # "We called the signals colocalized when (coloc H3+H4 ≥ 0.8 and H4∕H3 ≥ 2)" -Yi et al. (2019)
            report <- paste("Datasets colocalized")
        } else {
            report <- paste("Datasets NOT colocalized")
        }
    } else {
        report <- paste("Datasets NOT colocalized")
    }
    messager("+ COLOC::", report, "at: PP.H3 + PP.H4 >=", PP_threshold, " and PP.H3 / PP.H4 >= 2.")
    return(coloc_DT)
}
