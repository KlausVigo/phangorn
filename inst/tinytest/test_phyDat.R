# to test: phyDat.codon,
data(Laurasiatherian)
data(chloroplast)
set.seed(42)
tree <- rtree(10)
codon_align <- simSeq(tree, l=100, type = "CODON")


phy_matrix <- as.character(Laurasiatherian)
phy_df <- as.data.frame(Laurasiatherian)
phy_vec_dna <- phy_matrix[,1]
phy_vec_user <- sample(c("0","1"), 26, replace=TRUE)
names(phy_vec_user) <- letters
phy_dnabin <- as.DNAbin(Laurasiatherian)
phy_aabin <- as.AAbin(chloroplast)
phy_align <- phyDat2alignment(Laurasiatherian)

#test conversion work
expect_true(inherits(phy_matrix, "matrix"))
expect_true(inherits(phy_df, "data.frame"))
expect_true(inherits(phy_dnabin, "DNAbin"))
expect_true(inherits(phy_aabin, "AAbin"))
expect_true(inherits(phy_align, "alignment"))
expect_true(inherits(as.phyDat(phy_matrix), "phyDat"))
expect_true(inherits(as.phyDat(phy_df), "phyDat"))
expect_true(inherits(as.phyDat(phy_dnabin), "phyDat"))
# expect_equal(as.phyDat(phy_aabin), chloroplast))
expect_true(inherits(phyDat(phy_vec_dna), "phyDat"))
expect_true(inherits(phyDat(phy_vec_user, type="USER", levels = c("0","1")),
                    "phyDat"))
expect_true(inherits(as.phyDat(phy_dnabin), "phyDat"))
expect_true(inherits(as.phyDat(phy_align), "phyDat"))
expect_true(inherits(c2d <- codon2dna(codon_align), "phyDat"))
expect_equal(dna2codon(c2d), codon_align)


# test conversion with Biostrings
if(suppressPackageStartupMessages(requireNamespace('Biostrings'))){
  expect_true(inherits(MA_AA <- as.MultipleAlignment(chloroplast),
                  "AAMultipleAlignment"))
#  expect_equal(as.phyDat(MA_AA), chloroplast)
  expect_true(inherits(MA_DNA <- as.MultipleAlignment(Laurasiatherian),
                  "DNAMultipleAlignment"))
  expect_equal(as.phyDat(MA_DNA), Laurasiatherian)
}


# test subsetting and combining
expect_true(inherits(subset_1 <- subset(Laurasiatherian, select = 1:1000,
                                 site.pattern = FALSE), "phyDat"))
expect_true(inherits(subset_2 <- subset(Laurasiatherian, select = 1001:3179,
                                 site.pattern = FALSE), "phyDat"))
expect_true(inherits(lauraCbind1 <- cbind(subset_1, subset_2), "phyDat"))
expect_equal(baseFreq(lauraCbind1), baseFreq(Laurasiatherian))
expect_true(inherits(subset_3 <- subset(Laurasiatherian, select = 1:100),
                     "phyDat"))
expect_true(inherits(subset_4 <- subset(Laurasiatherian, select = 101:1605),
                  "phyDat"))
expect_equal(subset_1, Laurasiatherian[, 1:1000])
expect_error(subset(Laurasiatherian, 1:100), "subscript out of bounds")
expect_true(inherits(lauraCbind2 <- cbind(subset_3, subset_4), "phyDat"))
expect_equal(baseFreq(lauraCbind2), baseFreq(Laurasiatherian))


# test read and write
write.phyDat(Laurasiatherian, "tmp1.txt")
expect_true(inherits(laura <- read.phyDat("tmp1.txt"), "phyDat"))
expect_equal(laura, Laurasiatherian)
unlink("tmp1.txt")
write.phyDat(chloroplast, "tmp2.txt")
expect_true(inherits(chloro <- read.phyDat("tmp2.txt", type="AA"), "phyDat"))
# expect_equal(chloro, chloroplast) # changed to toupper
unlink("tmp2.txt")
write.phyDat(chloroplast, "tmp.fas", format="fasta")
expect_true(inherits(chloro_fas <- read.phyDat("tmp.fas", type="AA",
                                              format = "fasta"), "phyDat"))
# expect_equal(chloro_fas, chloroplast) # changed to toupper
unlink("tmp.fas")


# test removing duplicated sequences
tmp <- as.character(Laurasiatherian)
laura <- phyDat(rbind(phy_matrix, phy_matrix))
names(laura) <- paste0(names(laura), rep(c(1,2), each=47))

map1 <- map_duplicates(laura)
map2 <- map_duplicates(Laurasiatherian)
expect_null(map2)
expect_true(inherits(map1, "data.frame"))

