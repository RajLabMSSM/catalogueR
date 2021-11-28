# Interactive data table (when manually creating rmarkdown chunks)
createDT_html <- function(DF, caption = "", scrollY = 400) {
    htmltools::tagList(createDT(DF, caption, scrollY))
}
