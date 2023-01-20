data("Laurasiatherian")
set.seed(123)
trees_ultra <- bootstrap.phyDat(Laurasiatherian,
                                FUN=function(x)upgma(dist.ml(x)), bs=50)
trees_unrooted <- bootstrap.phyDat(Laurasiatherian,
                                   FUN=function(x)NJ(dist.ml(x)), bs=50)
tree_ultra <- allCompat(trees_ultra, rooted=TRUE)
tree_ultra_el <- add_edge_length(tree_ultra, trees_ultra, rooted=TRUE)
tree_unrooted <- allCompat(trees_unrooted, rooted=FALSE)
tree_unrooted_el <- add_edge_length(tree_unrooted, trees_unrooted, rooted=FALSE)

expect_null(tree_ultra$edge.length)
expect_equal(length(tree_ultra_el$edge.length), nrow(tree_ultra_el$edge))
expect_true(is.ultrametric(tree_ultra_el))
expect_null(tree_unrooted$edge.length)
expect_equal(length(tree_unrooted_el$edge.length), nrow(tree_unrooted_el$edge))
expect_false(is.ultrametric(tree_unrooted_el))
