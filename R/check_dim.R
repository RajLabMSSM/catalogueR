#' Check data dimensions
#' 
#' Check data dimensions and print to screen.
#' @keywords internal
check_dim <- function(df) {
    tryCatch({
        messager(
            "Data dimensions:",
            formatC(nrow(df), big.mark = ","), "rows",
            "x", formatC(ncol(df), big.mark = ","), "columns"
        )
    }, error = function(e){message(e);NULL})
}
