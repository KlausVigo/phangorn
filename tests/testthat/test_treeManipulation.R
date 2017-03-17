context("treeManipulation")

set.seed(42)
tree <- rtree(100, rooted=FALSE)
tree2 <- root(tree, 1, resolve.root = TRUE)


trees <- lapply(sample(10:500,50), function(x)tree <- rtree(x, rooted=FALSE) )

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
    ## ancestor, mrca, descendants
    expect_equal(mrca.phylo(tree, node=desc_108), 108L)
    kids_108 <- Descendants(tree, 108, "children")
    expect_equal(length(Descendants(tree, 101L, "all")), 197L)
    expect_equal(lengths(Descendants(tree2, 101L:199, "all")), 2 * lengths(prop.part(tree2)) - 2L)
    expect_equal(Ancestors(tree, kids_108, "parent"), rep(108L, length(kids_108)))
    expect_equal(Siblings(tree, kids_108[1], include.self=TRUE), kids_108)
})

test_that("allTrees", {
    ## allTrees
    expect_is(allTrees(6), "multiPhylo")
    expect_true(all(RF.dist(allTrees(6))>0))
})


test_that("nni", {
    ## nni
    expect_is(nni(tree), "multiPhylo")
    expect_true(all(RF.dist(nni(tree), tree)>0))
})


# TODO: check why rooted trees give error in development version
test_that("midpoint", {
    # topology stays the same
    expect_equal( max( sapply(trees, function(x)RF.dist(x,midpoint(x)))), 0) 
    # 2 * max(height) == max(cophenetic) 
    expect_equal( max( node.depth.edgelength(midpoint(tree)) *2) ,  max(cophenetic(tree)))              
})    
    


test_that("maxCladeCred", {
  tree <- rcoal(100)
  trees <- nni(tree)
  expect_equal(maxCladeCred(c(tree, trees)), tree)
})

