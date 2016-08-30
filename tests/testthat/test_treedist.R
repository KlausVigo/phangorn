context("treedist")

## generate data
set.seed(42)
tree <- rtree(10, FALSE)
trees <- nni(tree)
tree1 <- read.tree(text="(t5:1.0,(t4:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree2 <- read.tree(text="(t4:1.0,(t5:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree3 <- read.tree(text="(t5:1.0,t4:1.0,t3:1.0,(t1:1.0,t2:1.0):1.0);")

test_that("Robinson-Foulds distance", {
    ## check RF.dist and tree dist
    expect_that(RF.dist(tree1, tree1), is_equivalent_to(0))
    expect_that(RF.dist(tree1, tree2), is_equivalent_to(2))
    expect_that(RF.dist(tree1, tree3), is_equivalent_to(1))
    
    expect_true( all(RF.dist(tree, trees)==2) ) 
    expect_is(RF.dist(trees),"dist")
    
    expect_that(treedist(tree1, tree1)[1], is_equivalent_to(0))
    expect_that(treedist(tree1, tree2)[1], is_equivalent_to(2))
    expect_that(treedist(tree1, tree3)[1], is_equivalent_to(1))
})

test_that("Kuhner-Felsenstein distance (branch score difference)", {
    ## check RF.dist and tree dist
    expect_that(KF.dist(tree1, tree1), is_equivalent_to(0))
    expect_that(KF.dist(tree1, tree2), is_equivalent_to(sqrt(2)))
    expect_that(KF.dist(tree1, tree3), is_equivalent_to(1))
    
    expect_is(KF.dist(trees),"dist")
    
    expect_that(treedist(tree1, tree1)[2], is_equivalent_to(0))
    expect_that(treedist(tree1, tree2)[2], is_equivalent_to(sqrt(2)))
    expect_that(treedist(tree1, tree3)[2], is_equivalent_to(1))
    
    expect_equal(KF.dist(tree1, c(tree1, tree2, tree3)), c(0, sqrt(2), 1))
})

test_that("path distance", {
    ## check path dist
    expect_that(path.dist(tree1, tree1), is_equivalent_to(0))
    expect_equal(path.dist(trees)[1:13] , path.dist(trees[[1]], trees[2:14]))
    expect_is(path.dist(trees),"dist")
})


############################
# new tests from Michelle 
############################

# make simple trees with one pair of unmatched edges, unit branch lengths
tr1 <- read.tree(text="((A:1,B:1):1,C:1,(D:1,E:1):1);")
tr2 <- read.tree(text="((A:1,C:1):1,B:1,(D:1,E:1):1);")

# make a tree with same topology as tr1 but varied branch lengths
tr3 <- read.tree(text="((A:2,B:1):3,C:1,(D:1,E:2):1);")

test_that("Distance between known trees matches calculation by hand", {
    expect_equal(wRF.dist(tr1,tr2),2) # one pair of unmatched edges, unit branch lengths
    expect_equal(RF.dist(tr1,tr2),2) # one pair of unmatched edges (branch lengths irrelevant)
    expect_equal(wRF.dist(tr1,tr3),4) # same topology, different branch lengths
    expect_equal(RF.dist(tr1,tr3),0) # same topology (branch lengths irrelevant)
    expect_equal(wRF.dist(tr2,tr3),6) # one pair of unmatched edges, varied branch lengths
})

############################
# test that RF and wRF give same values for trees with every edge = 1
############################
test_that("When each tree has unit branch lengths, RF = wRF", {
    skip_on_cran()
    expect_equal(
        max(
            sapply(sample(10:500,50), function(x) { # test some random numbers of tips between 10 and 500
                trees <- rmtree(20, x, rooted=FALSE, br=1) # generate 20 unrooted trees with the given number of tips and unit branch lengths
                max(abs(RF.dist(trees) - wRF.dist(trees))) # find maximum abs difference between RF and wRF distance (expect 0)
            })), # find max of all these
        0) # expect equal to 0
    # now the same for comparison of one tree with many
    expect_equal(
        max(
            sapply(sample(10:500,50), function(x) { # test some random numbers of tips between 10 and 500
                trees <- rmtree(20, x, rooted=FALSE, br=1) # generate 20 unrooted trees with the given number of tips and unit branch lengths
                max(abs(RF.dist(trees[[1]], trees) - wRF.dist(trees[[1]], trees))) # find maximum abs difference between RF and wRF distance (expect 0)
            })), # find max of all these
        0) # expect equal to 0
})

#############################
# test sprdist from leomrtns 
#############################
test_that("SPR distance", {
    ## check spr dist
    set.seed(123)
    tree1 <- rtree(100, rooted = FALSE)
    tree2 <- rSPR(tree1, 1)
    trees <- rSPR(tree1, 1:5)
    expect_equal(sprdist(tree1, tree2)[[1]], 1)
    expect_equal(sprdist(tree1, tree2)[[3]], RF.dist(tree1, tree2))
    expect_equal(SPR.dist(tree1, trees), 1:5)
    expect_is(SPR.dist(trees), "dist")
})


