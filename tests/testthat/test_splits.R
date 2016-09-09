context("Test conversions")

## generate data
set.seed(1)
tree <- rtree(10, FALSE) 
tree2spl <- as.splits(tree)
spl2tree <- as.phylo(tree2spl)
dm <- cophenetic(tree2spl)
mat <- as.matrix(tree2spl)
Mat <- as.Matrix(tree2spl)
trees <- nni(tree)

test_that("splits ", {
    ## skip on CRAN
    skip_on_cran()
    
    ## check classes
    expect_is(as.splits(trees), "splits")
    expect_is(tree2spl, "splits")
    expect_is(spl2tree,"phylo")
    expect_is(dm,"dist")
    expect_is(mat, "matrix")
    expect_is(Mat, "Matrix")
    expect_equal(spl2tree , tree)
    
    spl <- allCircularSplits(6)
    spl <- phangorn:::oneWise(spl, 6)
    write.nexus.splits(spl, "tmp.nex")
    spl2 <- read.nexus.splits("tmp.nex")
    attr(spl2, "splitlabels") <- NULL
    attr(spl2, "weights") <- NULL
    class(spl2) <- "splits"
    expect_equal(spl2 , spl)
    unlink("tmp.nex")
})


test_that("networx ", {
     net1 <- neighborNet(dm)
     write.nexus.networx(net1, "tmp.nex")
     net2 <- read.nexus.networx("tmp.nex")
     net3 <- as.networx(tree)
     # delete some additional attributes
     net2$.plot <- net2$translate <- NULL
     attr(net1, "order") = NULL
     
     expect_is(net1, "networx")
     expect_is(net2, "networx")
     expect_is(net3, "networx")
#     expect_equal(net1, net2, tolerance=1e-6)
#     expect_equal(net3, net2, tolerance=1e-6)
     expect_equal(net1, net3)
     unlink("tmp.nex")
     cnet <- consensusNet(as.splits(trees))
     expect_is(cnet, "networx") 
     
     net1$edge.length <- cnet$edge.length <- cnet$edge.labels <- NULL
     attr(cnet, "order") = NULL
     expect_equal(cnet, net1)
})


