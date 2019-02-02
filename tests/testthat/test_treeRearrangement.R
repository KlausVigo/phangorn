context("tree_rearrangement")

set.seed(42)
tree_u <- rtree(100, rooted=FALSE)
tree_r <- rtree(100)


test_that("nni", {
    nni_trees_u <- nni(tree_u)
    nni_trees_r <- nni(tree_r)
    ## nni
    expect_is(nni(tree_u), "multiPhylo")
    expect_true(all(RF.dist(nni_trees_u, tree_u)>0))
    expect_true(length(nni_trees_u) == 194L)
    expect_true(length(nni_trees_r) == 196L)
    expect_true(all( RF.dist(tree_u, nni_trees_u) == 2))
    expect_true(median( RF.dist(tree_r, nni_trees_r) ) == 2)
})


r_nni <- rNNI(tree_u, 3, 100)


test_that("rNNI", {
    expect_true(length(r_nni) == 100L)
    expect_true(median( RF.dist(tree_u, r_nni) ) == 6)
    expect_true(median( RF.dist(r_nni[[1]], r_nni) ) == 12)
})

set.seed(42)
r_spr <- rSPR(tree_u, 3, 100)
test_that("allTrees", {
    expect_true(length(r_spr) == 100L)
    expect_true(median( SPR.dist(tree_u, r_spr) ) == 3)
    expect_true(median( SPR.dist(r_spr[[1]], r_spr) ) == 6)
})

