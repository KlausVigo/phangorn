context("speciesTree")

tr1 <- read.tree(text = "(((B:0.05,C:0.05):0.01,D:0.06):0.04,A:0.1);")
tr2 <- read.tree(text = "(((A:0.07,C:0.07):0.02,D:0.09):0.03,B:0.12);")
TR <- c(tr1, tr2)

start_tree <- read.tree(text = "(((B,C),D),A);")

x <- matrix(c("A", "B", "C", "D"), 4, 1)
rownames(x) <- c("A", "B", "C", "D") 
X <- phyDat(x, type="USER", levels = c("A", "B", "C", "D")) 


test_that("speciesTree", {
    ## check speciesTree
    st1 <- coalSpeciesTree(TR)
    st2 <- coalSpeciesTree(TR, X=X)
    st3 <- coalSpeciesTree(TR, sTree = start_tree)
    expect_equal(st1, st2)
    expect_equal(st1, st3)
})

