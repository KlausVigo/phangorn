fdir <- system.file("extdata/trees", package = "phangorn")
tmp <- read.csv(file.path(fdir,"H3N2_NA_20.csv"))
H3N2 <- read.phyDat(file.path(fdir,"H3N2_NA_20.fasta"), format="fasta")
dates <- setNames(tmp$numdate_given, tmp$name)

tree_unrooted <- candidate_tree(H3N2)
tree_ultra <- candidate_tree(H3N2, "ultra", tip.dates=dates)
tree_tipdated <- candidate_tree(H3N2, "tipdated", tip.dates=dates)

expect_false(is.rooted(tree_unrooted))
expect_true(is.rooted(tree_ultra))
expect_true(is.rooted(tree_tipdated))
expect_false(is.ultrametric(tree_tipdated))
expect_true(is.ultrametric(tree_ultra))

