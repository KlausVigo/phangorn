data("Laurasiatherian")
fun <- function(x) NJ(dist.hamming(x))
tree <- fun(Laurasiatherian)

bs_trees <- bootstrap.phyDat(Laurasiatherian, fun)
expect_inherits(bs_trees, "multiPhylo")

tree1 <- plotBS(tree, bs_trees, type="none")
tree2 <- plotBS(tree, bs_trees, type="none", method="TBE")

nnet <- neighborNet(dist.hamming(Laurasiatherian))

expect_true(inherits(tree1, "phylo"))
# transfer bootstrap should never be smaller than the standard one
expect_true(all(tree1$node.label[-1] <= tree2$node.label[-1]))


tree3 <- addConfidences(tree, bs_trees)
expect_equal(tree1, tree3)
nnet2 <- addConfidences(nnet, bs_trees)
expect_true(is.null(attr(nnet$splits, "confidences")))
expect_false(is.null(attr(nnet2$splits, "confidences")))

dat <- Laurasiatherian[sample(47, 10)] |> as.character() |> phyDat()
fit <- pml_bb(dat, "JC", rearrangement="NNI", control=pml.control(trace=0))
bs_trees_pml <- bootstrap.pml(fit, bs=10, rearrangement="NNI",
                              control=pml.control(trace=0))
expect_inherits(bs_trees_pml, "multiPhylo")


