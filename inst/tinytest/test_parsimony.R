data(yeast)
all_trees <- allTrees(8, tip.label = names(yeast))

tree1 <- read.tree(text = "((t1,t2),t3,t4);")
tree2 <- read.tree(text = "((t1,t3),t2,t4);")
trees <- .compressTipLabel(c(tree1, tree2))
dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER",
              levels=c("a","c","g","t"))

# test parsimony
expect_equal(fitch(tree1, dat), 1)
expect_equal(fitch(tree2, dat), 2)
expect_equal(fitch(trees, dat), c(1,2))
expect_equal(sankoff(tree1, dat), 1)
expect_equal(sankoff(tree2, dat), 2)
expect_equal(parsimony(tree1, dat), 1)


# test bab
#    all_trees <- allTrees(8, tip.label = names(yeast))
all_pars <- fitch(all_trees, yeast)
bab_tree <- bab(yeast, trace=0)
expect_equal(min(all_pars), fitch(bab_tree, yeast))


# test rearrangements
tree <- all_trees[[1]]
start <- fitch(tree, yeast)
bab_tree <- bab(yeast, trace=0)
best <- fitch(bab_tree, yeast)
best_fitch <- optim.parsimony(tree, yeast, rearrangements = "NNI", trace=0)
best_sankoff <- optim.parsimony(tree, yeast, method="sankoff",
                                rearrangements = "NNI", trace=0)
expect_equal(attr(best_fitch, "pscore"), attr(best_sankoff, "pscore"))


# test tree length
tree <- nj(dist.hamming(yeast))
pscore <- fitch(tree, yeast)
tree1 <- acctran(tree, yeast)
expect_equal(sum(tree1$edge.length), pscore)

tree2 <- rtree(100)
dat <- simSeq(tree2)
tree2 <- acctran(tree2, dat)
expect_equal(sum(tree2$edge.length), fitch(tree2,dat))



# test random.addition
ra_tree <- random.addition(yeast)
ratchet_tree <- pratchet(yeast, start=ra_tree, trace=0)
expect_true(attr(ra_tree, "pscore") >= attr(ratchet_tree, "pscore"))


