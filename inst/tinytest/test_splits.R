## generate data
set.seed(1)
tree <- rtree(10, FALSE)
tree2spl <- as.splits(tree)
pp <- prop.part(tree)
spl2tree <- as.phylo(tree2spl)
dm <- cophenetic(tree2spl)
mat <- as.matrix(tree2spl)
Mat <- as.Matrix(tree2spl)
trees <- nni(tree)

# test splits
## check classes
expect_true(inherits(as.splits(trees), "splits"))
expect_true(inherits(tree2spl, "splits"))
expect_true(inherits(spl2tree,"phylo"))
expect_true(inherits(dm,"dist"))
expect_true(inherits(mat, "matrix"))
expect_true(inherits(Mat, "Matrix"))
expect_equal(spl2tree , tree)
# test generics
c_spl <- c(tree2spl, tree2spl, tree2spl)
expect_equal(length(c_spl) , 3L*length(tree2spl))
expect_equal(length(unique(c_spl)) , length(tree2spl))
expect_equal(length(distinct.splits(c_spl)) , length(tree2spl))
spl <- allCircularSplits(6)
spl <- ONEwise(spl)
write.nexus.splits(spl, "tmp.nex")
spl2 <- read.nexus.splits("tmp.nex")
attr(spl2, "splitlabels") <- NULL
attr(spl2, "weights") <- NULL
class(spl2) <- "splits"
expect_equal(spl2 , spl)
# test conversion with prop.part
expect_true(inherits(as.splits(pp),"splits"))
expect_equal(pp, as.prop.part(as.splits(pp)))
expect_equivalent(as.splits(as.bitsplits(tree2spl)), tree2spl)
unlink("tmp.nex")


# test networx
net1 <- neighborNet(dm)
write.nexus.networx(net1, "tmp.nex")
net2 <- read.nexus.networx("tmp.nex")
net3 <- as.networx(tree)
# delete some additional attributes
net2$.plot <- net2$translate <- NULL
attr(net1, "order") <- NULL
expect_true(inherits(net1, "networx"))
expect_true(inherits(net2, "networx"))
expect_true(inherits(net3, "networx"))
expect_equal(net1, net2, use.edge.length = FALSE)
expect_equal(net3, net2, use.edge.length = FALSE)
expect_equal(net1, net3)
unlink("tmp.nex")

cnet <- consensusNet(trees)
expect_true(inherits(cnet, "networx"))
net1$edge.length <- cnet$edge.length <- cnet$edge.labels <- NULL
attr(cnet, "order") <- NULL
expect_equal(cnet, net1)
expect_equal(nrow(cnet$edge), length(as.splits(cnet)))




# test consensusNet
set.seed(1)
data("Laurasiatherian")
bs <- bootstrap.phyDat(Laurasiatherian,
                       FUN = function(x)nj(dist.hamming(x)), bs=50)
cnet <- consensusNet(bs, .2)
expect_true(inherits(cnet, "networx"))


spl <- allSplits(4)
net <- as.networx(spl)
expect_equal(Nnode(net), 8L)


