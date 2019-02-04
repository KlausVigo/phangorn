context("distances")

X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4, shape=0.5)

weights <- as.vector(1000*exp(fit$site))
attr(X, "weight") <- weights
dm <- cophenetic(tree)

Y <- phyDat(matrix(c("A", "C", "G", "T", "A", "C", "G", "A"), 2, 4,
                   dimnames=list(c("a", "b", NULL)), byrow=TRUE))

fun <- function(s) - 3/4 * log(1 - 4/3 * s)


test_that("dist.ml works properly", {
#    skip_on_cran()
    expect_that(dist.logDet(X), is_a("dist"))
    expect_that(dist.hamming(X), is_a("dist"))
    expect_that(dist.ml(X), is_a("dist"))
    expect_equal(as.matrix(dist.ml(X, k=4, shape=.5)), dm)
    expect_equal(as.matrix(dist.ml(Y)), as.matrix(fun(dist.hamming(Y))))
})


test_that("read/write of distances works", {
    skip_on_cran()
    writeDist(dm, "dm.txt")
    expect_equal(as.dist(dm), readDist("dm.txt"))
    unlink("dm.txt")
})



