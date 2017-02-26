library(testthat)

context("phangorn:::discrete.gamma")

test_that("149", {
    expected <- 1
    expect_equal({
        phangorn:::discrete.gamma(alpha = 1, k = 1L)
    }, expected)
})
