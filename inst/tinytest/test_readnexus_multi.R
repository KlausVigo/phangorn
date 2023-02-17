fcat <- function(..., file = zz) cat(..., file = file, sep = "",
                                     append = TRUE)
zz <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".nex")

tree <- rtree(10)
dat <- simSeq(tree, l=24)
write.phyDat(dat, file=zz, format="nexus")

fcat("BEGIN SETS;\n")
fcat("  Charset codon1 = 1-12/3;\n")
fcat("  Charset codon2 = 2-12/3;\n")
fcat("  Charset codon3 = 3-12/3;\n")
fcat("  Charset range = 16-18;\n")
fcat("  Charset range2 = 13-15 19-21;\n")
fcat("  Charset singles = 22 23 24;\n")
fcat("END;\n")

partitions <- setNames(c(4,4,4,3,6,3),
              c("codon1", "codon2", "codon3", "range", "range2", "singles"))

tmp <- read.nexus.partitions(zz)
length_part <- sapply(tmp, \(x)sum(attr(x, "weight")) )
expect_equal(partitions, length_part)


unlink(zz)
