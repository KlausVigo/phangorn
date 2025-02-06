
data("Laurasiatherian")

rooted_tree <- rcoal(47, names(Laurasiatherian))
unrooted_tree <- rtree(47, rooted=FALSE, names(Laurasiatherian))

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

#expect_equal, expect_equal_to_reference, expect_equivalent, expect_error,
#expect_false, expect_identical, expect_length, expect_match, expect_message,
#expect_null, expect_silent, expect_true, expect_warning, test_package
