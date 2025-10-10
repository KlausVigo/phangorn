trees <- rmtree(30, 5, rooted=FALSE)
trees <- .compressTipLabel(trees, ref=sort(trees$tip.label))

tree_u1 <- unique(trees)
tmp <- hash(trees)
tree_u2 <- trees[!duplicated(tmp)]

#expect_equal(attr(tree_u1, "old.index"),  match(tmp, unique(tmp)))
attr(tree_u1, "old.index") <- NULL
expect_equal(tree_u1, tree_u2)

#all.equal
# RF.dist  && dist.topo

