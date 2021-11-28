generate_token <- function(n = 1, len = 5) {
    a <- do.call(paste0, replicate(len, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}
