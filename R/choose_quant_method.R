# Try selecting the quant_method requested, but if not available select another
choose_quant_method <- function(di,
                                qm,
                                verbose = TRUE) {
    dataset_id <- quant_method <- NULL;
    meta <- eQTLcatalogue_list_datasets(verbose = FALSE)
    meta.sub <- data.frame(subset(meta, dataset_id == di))
    if (qm %in% unique(meta.sub$quant_method)) {
        meta.sub <- subset(meta.sub, quant_method == qm)
    } else {
        meta.sub <- meta.sub[1, ]
        messager("+ Selecting quant_method:", meta.sub$quant_method[1],
                 v = verbose)
    }
    return(meta.sub)
}
