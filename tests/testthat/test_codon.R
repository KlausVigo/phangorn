context("codon")


tree <- read.tree(text = "(((t1:.1,t2:.1):.1,t3:.2):.1,(t4:.1,t5:.1):.2);")
set.seed(42)
dat_dna <- simSeq(tree, l=300)
dat_codon <- dna2codon(dat_dna)

fit_F1x4 <- pml(tree, dat_codon, bf="F1x4")
fit_F3x4 <- pml(tree, dat_codon, bf="F3x4")
fit_GY <- pml(tree, dat_codon, bf="empirical")
