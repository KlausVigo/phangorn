context("treedist")

## generate data
set.seed(42)
tree <- rtree(10, FALSE)
trees <- nni(tree)
tree1 <- read.tree(text="(t5:1.0,(t4:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree2 <- read.tree(text="(t4:1.0,(t5:1.0,t3:1.0):1.0,(t1:1.0,t2:1.0):1.0);")
tree3 <- read.tree(text="(t5:1.0,t4:1.0,t3:1.0,(t1:1.0,t2:1.0):1.0);")

test_that("Robinson-Foulds distance", {
    ## skip_on_cran()
    
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
    ## skip_on_cran()
    
    ## check RF.dist and tree dist
    expect_that(KF.dist(tree1, tree1), is_equivalent_to(0))
    expect_that(KF.dist(tree1, tree2), is_equivalent_to(sqrt(2)))
    expect_that(KF.dist(tree1, tree3), is_equivalent_to(1))
    
    expect_is(KF.dist(trees),"dist")
    
    expect_that(treedist(tree1, tree1)[2], is_equivalent_to(0))
    expect_that(treedist(tree1, tree2)[2], is_equivalent_to(sqrt(2)))
    expect_that(treedist(tree1, tree3)[2], is_equivalent_to(1))
})

test_that("path distance", {
    ## skip_on_cran()
    
    ## check path dist
    expect_that(path.dist(tree1, tree1), is_equivalent_to(0))

    expect_is(path.dist(trees),"dist")
})

