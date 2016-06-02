context("treeManipulation")

set.seed(42)
tree <- rtree(100)

#Ancestors(x, node, type=c("all","parent"))
#Children(x, node)
#Siblings(x, node, include.self=FALSE)
#Descendants(x, node, type=c("tips","children","all"))
#Siblings(tree, 3)
#mrca.phylo(tree, 1:3)
#midpoint


desc_108 <- Descendants(tree, 108)[[1]]
node_108 <- mrca.phylo(tree, node=desc_108)

test_that("ancestor, mrca, descendants", {
    ## check RF.dist and tree dist
    expect_equal(mrca.phylo(tree, node=desc_108), 108)
})

test_that("allTrees, nni", {
    ## allTrees
    expect_is(allTrees(6), "multiPhylo")
    expect_true(all(RF.dist(allTrees(6))>0))
    ## nni
    expect_is(nni(tree), "multiPhylo")
    expect_true(all(RF.dist(nni(tree), tree)>0))
})





