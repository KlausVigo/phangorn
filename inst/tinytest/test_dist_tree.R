data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)

# check nnls functions
tree_nj <- NJ(dm)
tree_unj <- UNJ(dm)
tree_nnls_unj <- nnls.phylo(tree_unj, dm)
tree_nnls_nj <- nnls.phylo(tree_nj, dm)

tree_upgma <- upgma(dm)
tree_wpgma <- wpgma(dm)
tree_nnls_upgma <- nnls.phylo(tree_upgma, dm, method="ultrametric")
tree_nnls_wpgma <- nnls.phylo(tree_wpgma, dm, method="ultrametric")

expect_equal(tree_upgma, tree_nnls_upgma)
expect_false(all.equal(tree_wpgma, tree_nnls_wpgma))
expect_equal(tree_unj, tree_nnls_unj)
expect_false(all.equal(tree_nj, tree_nnls_nj))

# Test NNI
# tree_upgma_no_nni <- upgma(dm, NNI=FALSE)
# tree_upgma_nni <- upgma(dm, NNI=TRUE)

# expect_true(is.ultrametric(tree_upgma_nni))
# expect_true(sum(tree_upgma_no_nni$edge.length)>=sum(tree_upgma_nni$edge.length))

