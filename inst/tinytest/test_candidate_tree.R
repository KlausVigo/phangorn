# test candidate_tree
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

# test designTree and friends
dm <- dist.ml(H3N2)
tree_unj <- UNJ(dm)
tree_upgma <- upgma(dm)
tree_NNLS_unrooted <- nnls.tree(dm, tree_unj)
expect_equal(tree_NNLS_unrooted , tree_unj)
tree_NNLS_ultra <- nnls.tree(dm, tree_upgma, "ultra")
expect_equal(tree_NNLS_ultra , tree_upgma)

X_unj <- designTree(tree_unj)
lm_unrooted <- lm(dm ~ X_unj -1)

expect_equal(tree_unj$edge.length,lm_unrooted$coef, check.names = FALSE)

X_ultra <- designTree(tree_upgma, "ultra")
lm_ultra <- lm(dm ~ X_ultra -1)
coef_ultra <- numeric(37)
coef_ultra[20:37] <- lm_ultra$coef
el_ultra <- coef_ultra[tree_upgma$edge[,1]] - coef_ultra[tree_upgma$edge[,2]]
expect_equal(el_ultra, tree_upgma$edge.length)

