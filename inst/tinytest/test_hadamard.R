v <- 1:8
data(yeast)
dm <- dist.hamming(yeast)
# RY-coding
yeast_ry <- acgt2ry(yeast)
# delete ambiguous states
# dat4 <- phyDat(as.character(yeast), type="USER", levels=c("a","c", "g", "t"),
#                ambiguity=NULL)
# fit4 <- h4st(dat4)

# test Hadamard conjugation
# fast fft like multiplication
expect_true(inherits(H <- hadamard(3), "matrix"))
expect_equal(as.vector(H %*% v), fhm(v))
expect_true(inherits(spl_ry <- h2st(yeast_ry), "splits"))
expect_true(inherits(spl_dm <- distanceHadamard(dm), "splits"))

