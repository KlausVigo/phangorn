context("treeManipulation")


set.seed(42)
tree <- rtree(100)
trees_6 <- allTrees(6) 
trees_nni <- nni(tree)

#Ancestors(x, node, type=c("all","parent"))
#Children(x, node)
#Siblings(x, node, include.self=FALSE)
#Descendants(x, node, type=c("tips","children","all"))
#Ancestors(tree, 1:3, "all")
#Children(tree, 11)
#Descendants(tree, 11, "tips")
#Siblings(tree, 3)
#mrca.phylo(tree, 1:3)


desc_108 <- Descendants(tree, 108)[[1]]
node_108 <- mrca.phylo(tree, node=desc_108)

test_that("ancestor, mrca, descendants", {
    ## check RF.dist and tree dist
    expect_equal(mrca.phylo(tree, node=desc_108), 108)
})

test_that("allTrees, nni", {
    ## allTrees
    expect_is(trees_6, "multiPhylo")
    expect_true(all(RF.dist(trees_6)>0))
    ## nni
    expect_is(trees_nni, "multiPhylo")
    expect_true(all(RF.dist(trees_nni, tree)>0))
})





