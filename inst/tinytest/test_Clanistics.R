tree <- rtree(10)
x <- simSeq(tree, l = 5, type="USER", levels = c("red", "violet", "blue"))

# test that clanistics works properly
expect_inherits(getClans(tree), "matrix")
expect_inherits(getClips(tree), "matrix")
expect_inherits(getSlices(tree), "matrix")
expect_inherits(getDiversity(tree, x), "clanistics")
expect_inherits(getDiversity(tree, x), "data.frame")
