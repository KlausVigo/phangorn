library(testthat)

context("phangorn:::phyDat")

test_that("1", {
    expect_warning(expect_error({
        phangorn:::phyDat(data = structure(c("r", "a", "a", "a", 
            "a", "a", "y", "t", "t", "g", "g", "-", "g", "g", 
            "g", "a", "a", "a", "c", "t", "c", "-", "-", "c", 
            "c", "c", "c", "t", "t", "t", "c", "c", "?", "g", 
            "a", "g"), .Dim = c(3L, 12L), .Dimnames = list(c("t1", 
            "t2", "t3"), NULL)), levels = c("a", "c", "g", "t", 
            "-"), type = "USER", c("?", "n"))
    }, "argument is not interpretable as logical"), "the condition has length > 1 and only the first element will be used")
})
