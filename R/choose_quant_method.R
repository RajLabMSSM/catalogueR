# Try selecting the quant_method requested, but if not available select another
choose_quant_method <- function(ui,
                                qm,
                                verbose = TRUE) {
    meta <- eQTL_Catalogue.list_datasets(verbose = FALSE)
    meta.sub <- data.frame(subset(meta, unique_id == ui))
    if (qm %in% unique(meta.sub$quant_method)) {
        meta.sub <- subset(meta.sub, quant_method == qm)
    } else {
        meta.sub <- meta.sub[1, ]
        messager("+ Selecting quant_method:", meta.sub$quant_method[1], v = verbose)
    }
    return(meta.sub)
}
