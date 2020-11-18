tree <- rtree(10)
x <- simSeq(tree, l = 5, type="USER", levels = c("red", "violet", "blue"))

# test that clanistics works properly
expect_true(inherits(getClans(tree), "matrix"))
expect_true(inherits(getClips(tree), "matrix"))
expect_true(inherits(getSlices(tree), "matrix"))
expect_true(inherits(getDiversity(tree, x), "clanistics"))
expect_true(inherits(getDiversity(tree, x), "data.frame"))
