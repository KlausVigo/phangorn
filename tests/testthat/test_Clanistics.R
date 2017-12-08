context("Clanistics")

tree <- rtree(10)
x <- simSeq(tree, l=5, type="USER", levels = c("red", "violet", "blue"))

test_that("clanistics works properly", {
#    skip_on_cran()
    expect_is(getClans(tree), "matrix") 
    expect_is(getClips(tree), "matrix") 
    expect_is(getSlices(tree), "matrix")
    expect_is(getDiversity(tree, x), "clanistics")
    expect_is(getDiversity(tree, x), "data.frame")
})
