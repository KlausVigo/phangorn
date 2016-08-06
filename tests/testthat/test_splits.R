context("Test conversions")

## generate data
set.seed(1)
tree <- rtree(10, FALSE) 
tree2spl <- as.splits(tree)
spl2tree <- as.phylo(tree2spl)
dm <- cophenetic(tree2spl)
mat <- as.matrix(tree2spl)
Mat <- as.Matrix(tree2spl)


test_that("splits ", {
    ## skip on CRAN
    skip_on_cran()
    
    ## check classes
    expect_is(tree2spl, "splits")
    expect_is(spl2tree,"phylo")
    expect_is(dm,"dist")
    expect_is(mat, "matrix")
    expect_is(Mat, "Matrix")
    expect_equal(spl2tree , tree)
})


test_that("networx ", {
     net1 <- neighborNet(dm)
     write.nexus.networx(net1, "tmp.nex")
     net2 <- read.nexus.networx("tmp.nex")
     net3 <- as.networx(tree)
     # delete some additional attributes
     net2$.plot <- net2$translate <- NULL
     attr(net1, "order") = NULL
     expect_equal(net1, net2, tolerance=1e-6)
     expect_equal(net3, net2, tolerance=1e-6)
     expect_equal(net1, net3)
     unlink("tmp.nex")
})

