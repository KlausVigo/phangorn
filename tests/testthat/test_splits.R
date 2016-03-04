context("Test conversions")

#  methods(class="splits")
# addConfidences as.matrix      as.Matrix      as.networx     as.phylo       as.prop.part   cophenetic     c              print         [     
# as.splits.multiPhylo* as.splits.networx*    as.splits.phylo*      as.splits.prop.part* 

test_that("test: network <-> igraph", {
    ## skip on CRAN
    skip_on_cran()
    rm(list=ls())
    
    ## generate data
    set.seed(1)
    tree <- rtree(10, FALSE) 
    tree2spl <- as.splits(tree)
    spl2tree <- as.phylo(tree2spl)
    dm <- cophenetic(tree2spl)
    mat <- as.matrix(tree2spl)
    Mat <- as.Matrix(tree2spl)
    
    ## check classes
    expect_is(tree2spl, "splits")
    expect_is(spl2tree,"phylo")
    expect_is(dm,"dist")
    expect_is(mat, "matrix")
    expect_is(Mat, "Matrix")
    expect_equal(spl2tree , tree)

})
