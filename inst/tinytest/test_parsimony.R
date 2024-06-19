data(yeast)
all_trees <- allTrees(8, tip.label = names(yeast))

tree1 <- read.tree(text = "((t1,t2),t3,t4);")
tree2 <- read.tree(text = "((t1,t3),t2,t4);")
trees <- .compressTipLabel(c(tree1, tree2))
dat <- phyDat(c(t1="a", t2="a",t3="t",t4="t"), type="USER",
              levels=c("a","c","g","t"))

# TODO
# test sitewise pscores, ancestral and different states binary, DNA, AA
# sankoff + fitch
# CI, RI

# test parsimony
expect_equal(fitch(tree1, dat), 1)
expect_equal(fitch(tree2, dat), 2)
expect_equal(fitch(trees, dat), c(1,2))
expect_equal(sankoff(tree1, dat), 1)
expect_equal(sankoff(tree2, dat), 2)
expect_equal(parsimony(tree1, dat), 1)


# test bab
all_pars <- fitch(all_trees, yeast)
bab_tree <- bab(yeast, trace=0)
expect_equal(min(all_pars), fitch(bab_tree, yeast))


for(i in 1:10){
  tree100 <- rtree(100, rooted=FALSE)
  dat_2 <- simSeq(tree100, type="USER", levels=c("a", "b"))
  dat_3 <- simSeq(tree100, type="USER", levels=c("a", "b", "c"))
  dat_4 <- simSeq(tree100)
  pf_2 <- parsimony(tree100, dat_2, method = "fitch", site = "pscore")
  ps_2 <- parsimony(tree100, dat_2, method = "sankoff", site = "pscore")
  expect_equal(pf_2, ps_2)
  pf_3 <- parsimony(tree100, dat_3, method = "fitch", site = "pscore")
  ps_3 <- parsimony(tree100, dat_3, method = "sankoff", site = "pscore")
  expect_equal(pf_3, ps_3)
  pf_4 <- parsimony(tree100, dat_4, method = "fitch", site = "pscore")
  ps_4 <- parsimony(tree100, dat_4, method = "sankoff", site = "pscore")
  expect_equal(pf_4, ps_4)
  pvf_2 <- parsimony(tree100, dat_2, method = "fitch", site = "sitee")
  pvs_2 <- parsimony(tree100, dat_2, method = "sankoff", site = "site")
  expect_equal(pvf_2, pvs_2)
}




# test rearrangements
tree <- all_trees[[1]]
start <- fitch(tree, yeast)
bab_tree <- bab(yeast, trace=0)
best <- fitch(bab_tree, yeast)
best_fitch <- optim.parsimony(tree, yeast, rearrangements = "NNI", trace=0)
best_sankoff <- optim.parsimony(tree, yeast, method="sankoff",
                                rearrangements = "NNI", trace=0)
expect_equal(attr(best_fitch, "pscore"), attr(best_sankoff, "pscore"))


# test tree length
tree <- nj(dist.hamming(yeast))
pscore <- fitch(tree, yeast)
tree1 <- acctran(tree, yeast)
expect_equal(sum(tree1$edge.length), pscore)

tree2 <- rtree(100)
dat <- simSeq(tree2)
tree2 <- acctran(tree2, dat)
expect_equal(sum(tree2$edge.length), fitch(tree2,dat))


# test random.addition
ra_tree <- random.addition(yeast)
ratchet_tree <- pratchet(yeast, start=ra_tree, trace=0)
expect_true(attr(ra_tree, "pscore") >= attr(ratchet_tree, "pscore"))
trivial_tree <- pratchet(dat, trace=0, all=FALSE, minit = 10, maxit = 20)
expect_true(inherits(trivial_tree, "phylo"))


# test CI / RI
data(carnivora)
frm <- ~SuperFamily/Family/Genus/Species
tr <- as.phylo(frm, data = carnivora, collapse=FALSE)
tr$edge.length <- rep(1, nrow(tr$edge))

X <- matrix(0, 112, 7, dimnames = list(tr$tip.label, c("Canidae", "Felidae",
          "Ursidae", "Canidae_Felidae", "Canidae_Ursidae", "Felidae_Ursidae",
          "Canidae_Felidae_Ursidae")))
desc_canidae <- Descendants(tr, "Canidae")[[1]]
desc_felidae <- Descendants(tr, "Felidae")[[1]]
desc_ursidae <- Descendants(tr, "Ursidae")[[1]]

X[desc_canidae, c(1,4,5,7)] <- 1
X[desc_felidae, c(2,4,6,7)] <- 1
X[desc_ursidae, c(3,5,6,7)] <- 1

X <- phyDat(X, "USER", levels=c(0,1))

#col <- rep("black", 112)
#col[desc_felidae] <- "red"
#col[desc_canidae] <- "blue"
#col[desc_ursidae] <- "green"
#plot(tr, "f", tip.color=col, show.node=TRUE)

ci1 <- CI(tr, X, sitewise = TRUE)
ri1 <- RI(tr, X, sitewise = TRUE)
expect_true(all(ci1[1:3]==1))
expect_true(all(ri1[1:3]==1))

