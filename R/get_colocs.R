#' Run coloc on GWAS-QTL object
#'
#' @family coloc
#' @examples
#' library(dplyr)
#' data("BST1__Alasoo_2018.macrophage_IFNg")
#'
#' qtl.egene <- data.frame(BST1__Alasoo_2018.macrophage_IFNg)[, grep("*.QTL$|qtl_id|SNP", colnames(BST1__Alasoo_2018.macrophage_IFNg), value = TRUE)]
#' sorted_egenes <- qtl.egene %>%
#'     dplyr::group_by(gene.QTL) %>%
#'     dplyr::summarise(mean.P = mean(pvalue.QTL), min.P = min(pvalue.QTL)) %>%
#'     dplyr::arrange(min.P)
#' qtl.egene <- subset(qtl.egene, gene.QTL == sorted_egenes$gene.QTL[1])
#'
#' gwas.region <- data.frame(BST1__Alasoo_2018.macrophage_IFNg)[, grep("*.QTL$|qtl_id", colnames(BST1__Alasoo_2018.macrophage_IFNg), value = TRUE, invert = TRUE)]
#' coloc_res <- get_colocs(qtl.egene = qtl.egene, gwas.region = gwas.region)
get_colocs <- function(qtl.egene,
                       gwas.region,
                       merge_by_rsid = TRUE,
                       PP_threshold = .8,
                       method = "abf",
                       verbose = TRUE) {
    # http://htmlpreview.github.io/?https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/scripts/eQTL_API_usecase.html
    # Subset to overlapping SNPs only

    # Check for duplicated SNPs
    qtl.egene <- qtl.egene[!duplicated(qtl.egene$SNP), ]
    gwas.region <- gwas.region[!duplicated(gwas.region$SNP), ]
    if (all(c("N_cases", "N_controls") %in% colnames(gwas.region))) {
        message("Inferring `N` (effective sample size) from `N_cases` and `N_controls`")
        gwas.region$N <- get_sample_size(gwas.region, effective_ss = TRUE)$N
    }
    if (!"N" %in% colnames(gwas.region)) {
        stop("`N` column (effective sample size) was not detected in gwas.region. Required for coloc analysis.")
    }

    if (merge_by_rsid) {
        shared <- intersect(qtl.egene$SNP, gwas.region$SNP)
        eqtl_shared <- dplyr::filter(qtl.egene, SNP %in% shared) %>%
            dplyr::mutate(variant_id = SNP) %>%
            unique()
        gwas_shared <- dplyr::filter(gwas.region, SNP %in% shared) %>%
            dplyr::mutate(variant_id = SNP) %>%
            unique()
    } else {
        shared <- intersect(qtl.egene$position.QTL, gwas.region$POS)
        eqtl_shared <- dplyr::filter(qtl.egene, position.QTL %in% shared) %>%
            dplyr::mutate(variant_id = as.character(position.QTL))
        gwas_shared <- dplyr::filter(gwas.region, POS %in% shared) %>%
            dplyr::mutate(variant_id = as.character(POS))
    }

    #### Check MAF ####
    if (!"MAF" %in% colnames(gwas_shared)) {
        warning("`MAF` column not provided in GWAS data. Borrowing MAF from QTL data instead.")
        gwas_shared$MAF <- eqtl_shared$maf.QTL
    }

    if (length(shared) == 0) {
        if (verbose) {
            message("catalogueR:COLOC:: No SNPs shared between GWAS and QTL subsets.")
        }

        coloc_res <- list(
            summary = "No SNPs shared between GWAS and QTL subsets.",
            results = data.table::data.table(
                snp = NA,
                pvalues.df1 = NA,
                MAF.df1 = NA,
                V.df1 = NA,
                z.df1 = NA,
                r.df1 = NA,
                lABF.df1 = NA,
                V.df2 = NA,
                z.df2 = NA,
                r.df2 = NA,
                lABF.df2 = NA,
                internal.sum.lABF = NA
            ),
            Locus = gwas_shared$Locus[1]
        )
    } else {
        # RUN COLOC
        # QTL data
        eQTL_dataset <- list(
            snp = eqtl_shared$variant_id,
            pvalues = eqtl_shared$pvalue.QTL,
            # If log_OR column is full of NAs then use beta column instead
            beta = eqtl_shared$beta.QTL,
            N = (eqtl_shared$an.QTL)[1], # /2,
            MAF = eqtl_shared$maf.QTL,
            type = "quant"
        )
        if ("se.QTL" %in% colnames(eqtl_shared)) {
            eQTL_dataset$varbeta <- eqtl_shared$se.QTL^2
        }
        # GWAS data
        gwas_dataset <- list(
            snp = gwas_shared$variant_id,
            pvalues = gwas_shared$P,
            # If log_OR column is full of NAs then use beta column instead
            beta = gwas_shared$Effect,
            N = (gwas_shared$N)[1], # /2,
            MAF = gwas_shared$MAF,
            type = "cc",
            s = 0.5 # This is actually not used, because we already specified varbeta above.
        )
        if ("StdErr" %in% colnames(gwas_shared)) {
            gwas_dataset$varbeta <- gwas_shared$StdErr^2
        }

        # wrap <- ifelse(verbose, function(x)x, suppressMessages)

        # coloc::coloc.signals()
        if (method %in% c("single", "cond", "mask")) {
            # The updated coloc requires varbeta in both datasets
            color_res <- coloc::coloc.signals(
                dataset1 = eQTL_dataset,
                dataset2 = gwas_dataset,
                method = "mask",
                p12 = 1e-5
            )
        }
        if (tolower(method) == "abf") {
            coloc_res <- coloc::coloc.abf(
                dataset1 = eQTL_dataset,
                dataset2 = gwas_dataset,
                p12 = 1e-5
            ) # defaults
        }

        coloc_res$Locus <- gwas_shared$Locus[1]
        if (verbose) {
            report <- COLOC.report_summary(
                coloc.res = coloc_res,
                PP_threshold = PP_threshold
            )
        }
    }
    return(coloc_res)
}
