tree1 <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
tree2 <- read.tree(text = "((t1:0.3,t3:0.3):0.1,(t2:0.3,t4:0.3):0.1,t5:0.5);")

gene1 <- simSeq(tree1, l=200)
gene2 <- simSeq(tree1, l=200)
gene3 <- simSeq(tree1, l=200)
gene4 <- simSeq(tree2, l=200)
gene5 <- simSeq(tree2, l=200)

X <- cbind(gene1, gene2, gene3, gene4, gene5)
weight <- xtabs(~ index+genes,attr(X, "index"))

fit <- pml(tree1, X)
fit <- optim.pml(fit, control=pml.control(trace=0))

# test nni optimisation
sp <- pmlCluster( ~ edge + nni, fit, weight, p=1:3,
                     control=pml.control(epsilon=1e-08, maxit=10, trace=0))
#    expect_equal( sp$Partition, c(1,1,1,2,2))
expect_equal( all(sp$Partition ==  c(1,1,1,2,2)) ||
                      all(sp$Partition ==  c(2,2,2,1,1)) , TRUE)


