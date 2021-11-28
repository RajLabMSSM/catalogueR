#' Create an interactive data table with download buttons
#'
#' @param DF Data table.
#' @param caption Table caption.
#' @param scrollY Max height in pixels before scrolling.
#' @importFrom DT datatable
#' @export
createDT <- function(DF,
                     caption = "",
                     scrollY = 400) {
    data <- DT::datatable(DF,
        caption = caption, extensions = "Buttons",
        options = list(
            dom = "Bfrtip", buttons = c(
                "copy", "csv",
                "excel", "pdf", "print"
            ), scrollY = scrollY, scrollX = TRUE,
            scrollCollapse = TRUE, paging = FALSE, 
            columnDefs = list(list(
                className = "dt-center",
                targets = "_all"
            ))
        )
    )
    return(data)
}
