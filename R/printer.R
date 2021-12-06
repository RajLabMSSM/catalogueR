# Print
printer <- function(..., v = TRUE) {
    if (v) {
        msg <- paste(...)
        print(msg)
    }
}
