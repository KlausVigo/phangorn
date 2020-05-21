set.seed(42)
tree <- rtree(100, rooted=FALSE)
tree2 <- root(tree, 1, resolve.root = TRUE)


trees <- lapply(sample(10:500,50), function(x)tree <- rtree(x, rooted=FALSE) )


desc_108 <- Descendants(tree, 108)[[1]]
node_108 <- mrca.phylo(tree, node=desc_108)

# test ancestor, mrca, descendants
expect_equal(mrca.phylo(tree, node=desc_108), 108L)
expect_equal(mrca(tree), mrca.phylo(tree))
kids_108 <- Descendants(tree, 108, "children")
expect_equal(length(Descendants(tree, 101L, "all")), 197L)
expect_equal(lengths(Descendants(tree2, 101L:199, "all")),
                 2 * lengths(prop.part(tree2)) - 2L)
expect_equal(Ancestors(tree, kids_108, "parent"),
                 rep(108L, length(kids_108)))
expect_equal(Siblings(tree, kids_108[1], include.self=TRUE), kids_108)


# test allTrees
## allTrees
expect_true(inherits(allTrees(6), "multiPhylo"))
expect_true(all(RF.dist(allTrees(6))>0))


# TODO: check why rooted trees give error in development version
# test midpoint
# topology stays the same
expect_equal( max( sapply(trees, function(x)RF.dist(x,midpoint(x)))), 0)
# 2 * max(height) == max(cophenetic)
expect_equal( max( node.depth.edgelength(midpoint(tree)) *2),
                  max(cophenetic(tree)))


# test maxCladeCred
tree <- rcoal(100)
trees <- nni(tree)
expect_equal(maxCladeCred(c(tree, trees)), tree)
tree <- rtree(100, rooted = FALSE)
trees <- nni(tree)
expect_equal(tree, allCompat(trees), use.edge.length = FALSE)


# test add.tips
tree <- rcoal(20)
tree$tip.label <- paste0("t", 1:Ntip(tree))
tree$node.label <- paste0("n", 1:Nnode(tree))
tree1 <- add.tips(tree, c("A", "B", "C", "D"), c("t5", "t10", "n5", "n10"))
tree2 <- add.tips(tree, c("A", "B", "C", "D"), c(5, 10, 25, 30))
expect_equal(RF.dist(tree1, tree2), 0)
expect_false(is.binary(tree1))
expect_false(is.binary(tree2))


# test plotBS
set.seed(1)
data("Laurasiatherian")
bs <- bootstrap.phyDat(Laurasiatherian,
                       FUN = function(x)NJ(dist.hamming(x)), bs=50)
tree <- NJ(dist.hamming(Laurasiatherian))
treeBS <- plotBS(tree, bs, type="none")
expect_equal(length(treeBS$node.label), treeBS$Nnode)
expect_true(all(treeBS$node.label >=0) )
expect_true(all(treeBS$node.label <= 100) )
