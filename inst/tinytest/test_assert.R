data("Laurasiatherian")

set.seed(42)

rooted_tree <- rcoal(47, names(Laurasiatherian))
unrooted_tree <- rtree(47, rooted = FALSE, names(Laurasiatherian))
trees <- rmtree(10,10)

fit <- pml(unrooted_tree, Laurasiatherian)


tree_2 <- di2multi(rooted_tree, tol = 1e-02)
trees_2 <- c(tree_2, rooted_tree, unrooted_tree)



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

expect_silent(assert_multiPhylo(trees))
expect_error(assert_multiPhylo(rooted_tree))

expect_silent(assert_treeish(unrooted_tree))
expect_silent(assert_treeish(trees))
expect_error(assert_treeish(NULL))
expect_silent(assert_treeish(NULL, null.ok = TRUE))


# test clean_phylo
expect_false(is.binary(tree_2))
expect_false(is.binary(trees_2) |> all())
expect_true(clean_phylo(tree_2, multi2di = TRUE) |> is.binary())
expect_true(clean_multiPhylo(trees_2, multi2di = TRUE) |> is.binary() |> all())
