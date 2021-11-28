# Print
messager <- function(..., v = TRUE) {
    if (v) {
        print(paste(...))
    }
}
