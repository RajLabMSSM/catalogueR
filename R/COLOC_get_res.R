#' Run coloc on GWAS-QTL object
#'
#' Run coloc on GWAS-QTL object.
#' @source
#' \code{
#' library(dplyr)
#' paths <- catalogueR::eQTLcatalogue_example_queries()
#' query <- paths["BST1__Alasoo_2018.macrophage_IFNg"]
#'
#' qtl.egene <- data.frame(
#'     query
#'     )[, grep("*.QTL$|qtl_id|SNP",
#'              colnames(query), value = TRUE)]
#' sorted_egenes <- qtl.egene |>
#'     dplyr::group_by(gene.QTL) |>
#'     dplyr::summarise(mean.P = mean(pvalue.QTL), min.P = min(pvalue.QTL)) |>
#'     dplyr::arrange(min.P)
#' qtl.egene <- subset(qtl.egene, gene.QTL == sorted_egenes$gene.QTL[1])
#'
#' gwas.region <- data.frame(
#'     query
#'     )[, grep("*.QTL$|qtl_id",
#'              colnames(query),
#'              value = TRUE, invert = TRUE)]
#' #### Run ####
#' coloc_res <- catalogueR:::COLOC_get_res(qtl.egene = qtl.egene,
#'                                      gwas.region = gwas.region)
#' }
#' @family coloc
#' @keywords internal
#' @importFrom dplyr filter
#' @importFrom coloc coloc.signals coloc.abf
COLOC_get_res <- function(qtl.egene,
                       gwas.region,
                       merge_by_rsid = TRUE,
                       coloc_thresh = .8,
                       method = "abf",
                       verbose = TRUE) {
  
    SNP <- position.QTL <- POS <- NULL;
    method <- tolower(method)[[1]]
    #### Check for duplicated SNPs ####
    qtl.egene <- qtl.egene[!duplicated(qtl.egene$SNP), ]
    gwas.region <- gwas.region[!duplicated(gwas.region$SNP), ]
    if (!"N" %in% colnames(gwas.region)) {
        stp <- paste(
            "`N` column (effective sample size) was not detected",
            "in gwas.region. Required for coloc analysis."
        )
        stop(stp)
    }
    #### Merge by SNP IDs ####
    if (merge_by_rsid) {
        shared <- intersect(qtl.egene$SNP, gwas.region$SNP)
        eqtl_shared <- dplyr::filter(qtl.egene, SNP %in% shared) |>
            dplyr::mutate(variant_id = SNP) |>
            unique()
        gwas_shared <- dplyr::filter(gwas.region, SNP %in% shared) |>
            dplyr::mutate(variant_id = SNP) |>
            unique()
    } else {
        shared <- intersect(qtl.egene$position.QTL, gwas.region$POS)
        eqtl_shared <- dplyr::filter(qtl.egene, position.QTL %in% shared) |>
            dplyr::mutate(variant_id = as.character(position.QTL))
        gwas_shared <- dplyr::filter(gwas.region, POS %in% shared) |>
            dplyr::mutate(variant_id = as.character(POS))
    }
    #### Check MAF ####
    check_maf_out <- check_maf(
        gwas_shared = gwas_shared,
        eqtl_shared = eqtl_shared,
        verbose = verbose
    )
    gwas_shared <- check_maf_out$gwas_shared
    eqtl_shared <- check_maf_out$eqtl_shared

    if (length(shared) == 0) {
        messager("No SNPs shared between GWAS and QTL subsets.", v = verbose)
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
        #### Run coloc ####
        # QTL data
        eQTL_dataset <- list(
            snp = eqtl_shared$variant_id,
            pvalues = eqtl_shared$pvalue.QTL,
            # If log_OR column is full of NAs then use beta column instead
            beta = eqtl_shared$beta.QTL,
            N = (eqtl_shared$an.QTL)[1], # /2,
            MAF = as.numeric(eqtl_shared$maf.QTL),
            type = "quant"
        )
        if ("se.QTL" %in% colnames(eqtl_shared)) {
            eQTL_dataset$varbeta <- eqtl_shared$se.QTL^2
        }
        # GWAS data
        query_datset <- list(
            snp = gwas_shared$variant_id,
            pvalues = gwas_shared$P,
            # If log_OR column is full of NAs then use beta column instead
            beta = gwas_shared$Effect,
            N = (gwas_shared$N)[1], # /2,
            MAF = as.numeric(gwas_shared$MAF),
            type = "cc",
            # `s=` is actually not used, 
            # because we already specified varbeta above.
            s = 0.5 
        )
        if ("StdErr" %in% colnames(gwas_shared)) {
            query_datset$varbeta <- gwas_shared$StdErr^2
        } 
        if (method %in% c("single", "cond", "mask")) {
            # The updated coloc requires varbeta in both datasets
            coloc_res <- coloc::coloc.signals(
                dataset1 = eQTL_dataset,
                dataset2 = query_datset,
                method = method,
                MAF = NULL
            )
        }
        if (method == "abf") { 
                coloc_res <- suppressWarnings(
                  coloc::coloc.abf(
                    dataset1 = eQTL_dataset,
                    dataset2 = query_datset,
                    MAF = NULL)
                  ) 
        }

        coloc_res$Locus <- gwas_shared$Locus[1]
        if (isTRUE(verbose)) {
            report <- COLOC.report_summary(
                coloc.res = coloc_res,
                coloc_thresh = coloc_thresh
            )
        }
    }
    return(coloc_res)
}

#### Deprecation function #####
get_colocs <- function(...){
  .Deprecated("COLOC_get_res")
  COLOC_get_res(...)
}
