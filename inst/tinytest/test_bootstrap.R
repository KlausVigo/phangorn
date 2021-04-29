data("Laurasiatherian")
fun <- function(x) NJ(dist.hamming(x))
tree <- fun(Laurasiatherian)

bs_trees <- bootstrap.phyDat(Laurasiatherian, fun)
tree1 <- plotBS(tree, bs_trees, type="none")
tree2 <- plotBS(tree, bs_trees, type="none", method="TBE")


expect_true(inherits(tree1, "phylo"))
# transfer bootstrap should never be smaller than the standard one
expect_true(all(tree1$node.label[-1] <= tree2$node.label[-1]))

