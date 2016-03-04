context("treedist")

test_that("Robinson-Foulds distance", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())
    
    ## generate data
    tree <- rtree(10, FALSE)
    trees <- nni(tree)
    tree1 = read.tree(text="(t5:1.0,(t4:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
    tree2 = read.tree(text="(t4:1.0,(t5:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
    tree3 = read.tree(text="(t5:1.0,t4:1.0,t3:1.0,(t1:1.0,t2:1.0):1.0);")
    
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