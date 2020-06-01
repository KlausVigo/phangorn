X <- allSitePattern(5)
tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
fit <- pml(tree,X, k=4, shape=0.5)

weights <- as.vector(1000*exp(fit$siteLik))
attr(X, "weight") <- weights
dm <- cophenetic(tree)

Y <- phyDat(matrix(c("A", "C", "G", "T", "A", "C", "G", "A"), 2, 4,
                   dimnames=list(c("a", "b", NULL)), byrow=TRUE))

fun <- function(s) - 3/4 * log(1 - 4/3 * s)

data(woodmouse)

# test dist.ml
    expect_true(inherits(dist.logDet(X), "dist"))
    expect_true(inherits(dist.hamming(X), "dist"))
    expect_true(inherits(dist.ml(X), "dist"))
    expect_equal(as.matrix(dist.ml(X, k=4, shape=.5)), dm)
    expect_equal(as.matrix(dist.ml(Y)), as.matrix(fun(dist.hamming(Y))))
    expect_equivalent(dist.dna(woodmouse, "JC", pairwise.deletion = FALSE),
                      dist.ml(woodmouse, exclude = "all"))
    expect_equivalent(dist.dna(woodmouse, "JC", pairwise.deletion = TRUE),
                      dist.ml(woodmouse, exclude = "pairwise"))
    expect_equivalent( dist.dna(woodmouse, "raw"),
                       dist.hamming(woodmouse, exclude="all"))
    expect_equivalent( dist.dna(woodmouse, "raw", pairwise.deletion = TRUE),
                       dist.hamming(woodmouse, exclude="pairwise"))
    expect_equivalent( dist.dna(woodmouse, "N"),
                       dist.hamming(woodmouse, exclude="all", ratio = FALSE))
    expect_equivalent( dist.dna(woodmouse, "N", pairwise.deletion = TRUE),
                       dist.hamming(woodmouse, exclude="pairwise",
                                    ratio = FALSE))


# test read/write of distances
    # phylip
    dm <- as.dist(dm)
    writeDist(dm, "dm.txt")
    expect_equal(dm, readDist("dm.txt"))
    #nexus
    writeDist(dm, "dm.txt", format="nexus", upper=TRUE)
    expect_equal(dm, readDist("dm.txt", format="nexus"))
    writeDist(dm, "dm.txt", format="nexus", upper=FALSE)
    expect_equal(dm, readDist("dm.txt", format="nexus"))
    unlink("dm.txt")




