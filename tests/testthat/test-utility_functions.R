test_that("messager prints when v=TRUE", {
    expect_message(
        catalogueR:::messager("hello", "world", v = TRUE),
        "hello world"
    )
})

test_that("messager is silent when v=FALSE", {
    expect_no_message(
        catalogueR:::messager("hello", v = FALSE)
    )
})

test_that("fix_ftp converts ftp to http", {
    expect_equal(
        catalogueR:::fix_ftp("ftp://example.com/file.gz"),
        "http://example.com/file.gz"
    )
    # Already http - should remain unchanged
    expect_equal(
        catalogueR:::fix_ftp("http://example.com/file.gz"),
        "http://example.com/file.gz"
    )
    # https should remain unchanged
    expect_equal(
        catalogueR:::fix_ftp("https://example.com/file.gz"),
        "https://example.com/file.gz"
    )
})

test_that("timeit returns function result", {
    result <- catalogueR:::timeit(42 + 1)
    expect_equal(result, 43)

    result2 <- catalogueR:::timeit({
        x <- 10
        x * 2
    })
    expect_equal(result2, 20)
})

test_that("message_parallel produces output", {
    # message_parallel uses system echo, just verify it runs without error
    expect_no_error(
        catalogueR:::message_parallel("test message")
    )
})
