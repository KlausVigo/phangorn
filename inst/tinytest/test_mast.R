context("mast")

## generate data
set.seed(42)
tree1 <- rtree(100)
tree2 <- rSPR(tree1, 5)
tips <- mast(tree1, tree2, tree = FALSE)
mast_tree <- mast(tree1, tree2)

tip_label <- tree1$tip.label

tips_to_delete <- setdiff(tip_label, tips)

tree1_drop <- drop.tip(tree1, tips_to_delete)
tree2_drop <- drop.tip(tree2, tips_to_delete)


test_that("maximum agreement subtree (MAST)", {
    ## common subtrees should be identical
    expect_equal(RF.dist(tree1_drop, tree2_drop), 0)
    expect_equal(RF.dist(tree1_drop, mast_tree), 0)
    expect_equal(RF.dist(tree2_drop, mast_tree), 0)
})
