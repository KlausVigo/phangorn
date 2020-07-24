## generate data
library(ape)
set.seed(42)
tree <- rtree(10, FALSE)
trees <- nni(tree)
tree1 <- read.tree(text="(t5:1.0,(t4:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree2 <- read.tree(text="(t4:1.0,(t5:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree3 <- read.tree(text="(t5:1.0,t4:1.0,t3:1.0,(t1:1.0,t2:1.0):1.0);")

# Robinson-Foulds distance
## check RF.dist and tree dist
expect_equivalent(RF.dist(tree1, tree1), 0)
expect_equivalent(RF.dist(tree1, tree2), 2)
expect_equivalent(RF.dist(tree1, tree3), 1)

expect_true( all(RF.dist(tree, trees)==2) )
expect_true(inherits(RF.dist(trees),"dist"))

expect_equivalent(treedist(tree1, tree1)[1], 0)
expect_equivalent(treedist(tree1, tree2)[1], 2)
expect_equivalent(treedist(tree1, tree3)[1], 1)


# Kuhner-Felsenstein distance (branch score difference)
## check RF.dist and tree dist
expect_equivalent(KF.dist(tree1, tree1), 0)
expect_equivalent(KF.dist(tree1, tree2), sqrt(2))
expect_equivalent(KF.dist(tree1, tree3), 1)

expect_true(inherits(KF.dist(trees),"dist"))

expect_equivalent(treedist(tree1, tree1)[2], 0)
expect_equivalent(treedist(tree1, tree2)[2], sqrt(2))
expect_equivalent(treedist(tree1, tree3)[2], 1)

expect_equal(KF.dist(tree1, c(tree1, tree2, tree3)), c(0, sqrt(2), 1))


# path distance
## check path dist
expect_equivalent(path.dist(tree1, tree1), 0)
expect_equal(path.dist(trees)[1:13] , path.dist(trees[[1]], trees[2:14]))
expect_true(inherits(path.dist(trees),"dist"))



############################
# new tests from Michelle
############################

# make simple trees with one pair of unmatched edges, unit branch lengths
tr1 <- read.tree(text="((A:1,B:1):1,C:1,(D:1,E:1):1);")
tr2 <- read.tree(text="((A:1,C:1):1,B:1,(D:1,E:1):1);")

# make a tree with same topology as tr1 but varied branch lengths
tr3 <- read.tree(text="((A:2,B:1):3,C:1,(D:1,E:2):1);")

# Distance between known trees matches calculation by hand
# one pair of unmatched edges, unit branch lengths
expect_equal(wRF.dist(tr1,tr2),2)
# one pair of unmatched edges (branch lengths irrelevant)
expect_equal(RF.dist(tr1,tr2),2)
# same topology, different branch lengths
expect_equal(wRF.dist(tr1,tr3),4)
# same topology (branch lengths irrelevant)
expect_equal(RF.dist(tr1,tr3),0)
# one pair of unmatched edges, varied branch lengths
expect_equal(wRF.dist(tr2,tr3),6)


############################
# test that RF and wRF give same values for trees with every edge = 1
############################
# When each tree has unit branch lengths, RF = wRF", {
  expect_equal(
    max(
      # test some random numbers of tips between 10 and 500
      sapply(sample(10:500,50), function(x) {
      # generate 20 unrooted trees with the given number of tips and unit branch
      # lengths
      trees <- rmtree(20, x, rooted=FALSE, br=1)
      # find maximum abs difference between RF and wRF distance (expect 0)
      max(abs(RF.dist(trees) - wRF.dist(trees)))
      })), # find max of all these
      0) # expect equal to 0
    # now the same for comparison of one tree with many
  expect_equal(
    max(
      # test some random numbers of tips between 10 and 500
      sapply(sample(10:500,50), function(x) {
      # generate 20 unrooted trees with the given number of tips and unit branch
      # lengths
      trees <- rmtree(20, x, rooted=FALSE, br=1)
      # find maximum abs difference between RF and wRF distance (expect 0)
      max(abs(RF.dist(trees[[1]], trees) - wRF.dist(trees[[1]], trees)))
      })), # find max of all these
      0) # expect equal to 0


#############################
# test sprdist from leomrtns
#############################
# SPR distance
set.seed(123)
tree1 <- rtree(100, rooted = FALSE)
tree2 <- rSPR(tree1, 1)
trees <- rSPR(tree1, 1:5)
expect_equal(sprdist(tree1, tree2)[[1]], 1)
expect_equal(sprdist(tree1, tree2)[[3]], RF.dist(tree1, tree2))
expect_equal(SPR.dist(tree1, trees), 1:5)
expect_true(inherits(SPR.dist(trees), "dist"))

# Fix to issue #97 on github
tr1 <- structure(list(edge = structure(c(11L, 11L, 10L, 10L, 9L, 9L, 8L, 8L, 7L,
                                        7L, 2L, 6L, 5L, 11L, 4L, 10L, 3L, 9L,
                                        1L, 8L), .Dim = c(10L, 2L)),
                     tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"),
                     Nnode = 5), class = "phylo")

tr3 <- structure(list(edge = structure(c(9L, 9L, 11L, 11L, 10L, 10L, 8L, 8L, 7L,
                                         7L, 1L, 2L, 4L, 5L, 6L, 11L, 3L, 9L,
                                         8L, 10L), .Dim = c(10L, 2L)),
                      tip.label = c("t1", "t2", "t3", "t4", "t5", "t6"),
                      Nnode = 5L), class = "phylo")
expect_equal(SPR.dist(tr1, tr3), SPR.dist(tr3, tr1))
