
data("Laurasiatherian")

rooted_tree <- rcoal(47, names(Laurasiatherian))
unrooted_tree <- rtree(47, rooted=FALSE, names(Laurasiatherian))
trees <- rmtree(10,10)

fit <- pml(unrooted_tree, Laurasiatherian)

expect_error( assert_phylo(Laurasiatherian) )
expect_silent( assert_phylo(rooted_tree) )
expect_silent( assert_phylo(rooted_tree, is_rooted = TRUE) )
expect_silent( assert_phylo(rooted_tree, is_ultrametric = TRUE) )
expect_silent( assert_phylo(unrooted_tree) )
expect_error( assert_error(unrooted_tree, is_rooted = TRUE) )
expect_error( assert_error(unrooted_tree, is_ultrametric = TRUE) )

expect_silent( assert_phyDat(Laurasiatherian) )
expect_error( assert_phyDat(rooted_tree) )

expect_silent( assert_pml(fit) )
expect_error( assert_pml(rooted_tree) )


expert_silent(assert_treeish(unrooted_tree))
expert_silent(assert_treeish(trees))
expert_error(assert_treeish(NULL))
expert_silent(assert_treeish(NULL, null.ok=TRUE))
