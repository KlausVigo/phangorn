tree <- rtree(50, rooted=FALSE)

trees_simple <- nni(tree)

trees <- rNNI(tree, sample(10, 100, replace = TRUE))
trees <- .uncompressTipLabel(trees)
labels <- paste0("t", 1:50)
trees2 <- lapply(trees, function(x)drop.tip(x, sample(labels, 10)))
class(trees2) <- "multiPhylo"

# superTree
simple_superTree <- superTree(trees_simple, rooted=FALSE)
difficult_superTree <- superTree(trees2, rooted=FALSE)
expect_equal(RF.dist(simple_superTree, tree) , 0)
expect_equal(RF.dist(difficult_superTree, tree) , 0)
rf_superTree <- superTree(trees, method="RF")
spr_superTree <- superTree(trees, method="SPR")
expect_true(attr(rf_superTree, "score") > attr(spr_superTree, "score"))




