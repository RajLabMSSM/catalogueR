check_dim <- function(df) {
    try({
        messager("Data dimensions:", nrow(df), "x", ncol(df))
    })
}
