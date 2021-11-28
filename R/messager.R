messager <- function(..., v = T) {
    if (v) {
        message(paste(...))
    }
}
