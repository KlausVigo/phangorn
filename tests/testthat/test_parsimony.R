context("parsimony")

data(yeast)
all_trees <- allTrees(8, tip.label = names(yeast))

tree1 <- read.tree(text = "((t1,t2),t3,t4);")
tree2 <- read.tree(text = "((t1,t3),t2,t4);")
trees <- .compressTipLabel(c(tree1, tree2))
dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER", 
              levels=c("a","c","g","t"))

test_that("parsimony works properly", {
##    skip_on_cran()
    expect_that(fitch(tree1, dat), equals(1))
    expect_that(fitch(tree2, dat), equals(2))
    expect_that(fitch(trees, dat), equals(c(1,2)))
    expect_that(sankoff(tree1, dat), equals(1))
    expect_that(sankoff(tree2, dat), equals(2))
    expect_that(parsimony(tree1, dat), equals(1))
})

test_that("bab works properly", {  
    skip_on_cran()
#    all_trees <- allTrees(8, tip.label = names(yeast))
    all_pars <- fitch(all_trees, yeast)
    bab_tree <- bab(yeast, trace=0)
    expect_equal(min(all_pars), fitch(bab_tree, yeast))
})

test_that("rearrangements works properly", {  
    skip_on_cran()
    tree <- all_trees[[1]]
    start <- fitch(tree, yeast)
    bab_tree <- bab(yeast, trace=0)
    best <- fitch(bab_tree, yeast)
    best_fitch <- optim.parsimony(tree, yeast, rearrangements = "NNI", trace=0)
    best_sankoff <- optim.parsimony(tree, yeast, method="sankoff", 
                                    rearrangements = "NNI", trace=0)
    expect_equal(attr(best_fitch, "pscore"), attr(best_sankoff, "pscore"))
})


test_that("tree length works properly", {  
    skip_on_cran()
    tree <- nj(dist.hamming(yeast))
    pscore <- fitch(tree, yeast)
    tree1 <- acctran(tree, yeast)
    expect_equal(sum(tree1$edge.length), pscore)
    tree2 <- rtree(100)
    dat <- simSeq(tree2)
    tree2 <- acctran(tree2, dat)
    expect_equal(sum(tree2$edge.length), fitch(tree2,dat))
})



