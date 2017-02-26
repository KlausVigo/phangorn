library(testthat)

context("phangorn:::fast.table")

test_that("0", {
    expected <- structure(list(index = c(1, 2, 3, 4, 5, 2, 6, 
        7, 8, 9, 10, 11), weights = c(1L, 2L, 1L, 1L, 1L, 1L, 
        1L, 1L, 1L, 1L, 1L), data = structure(list(t1 = c("r", 
        "a", "y", "g", "g", "c", "-", "c", "t", "c", "g"), t2 = c("a", 
        "a", "t", "g", "g", "t", "-", "c", "t", "c", "a"), t3 = c("a", 
        "a", "t", "-", "g", "c", "c", "c", "t", "?", "g")), .Names = c("t1", 
        "t2", "t3"), row.names = c(1L, 2L, 3L, 4L, 5L, 7L, 8L, 
        9L, 10L, 11L, 12L), class = "data.frame")), .Names = c("index", 
        "weights", "data"))
    expect_equal({
        phangorn:::fast.table(data = structure(list(t1 = c("r", 
            "a", "y", "g", "g", "a", "c", "-", "c", "t", "c", 
            "g"), t2 = c("a", "a", "t", "g", "g", "a", "t", "-", 
            "c", "t", "c", "a"), t3 = c("a", "a", "t", "-", "g", 
            "a", "c", "c", "c", "t", "?", "g")), .Names = c("t1", 
            "t2", "t3"), row.names = c(NA, -12L), class = "data.frame"))
    }, expected)
})
