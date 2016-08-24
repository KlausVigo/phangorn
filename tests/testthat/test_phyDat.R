context("conversion_and_subsetting")
# to test: phyDat.codon, 
data(Laurasiatherian)
data(chloroplast)
set.seed(42)
tree <- rtree(10)
codon_align <- simSeq(tree, l=100, type = "CODON")


phy_matrix <- as.character(Laurasiatherian)
phy_df <- as.data.frame(Laurasiatherian)
phy_vec <- phy_matrix[,1]
phy_dnabin <- as.DNAbin(Laurasiatherian)
phy_align <- phyDat2alignment(Laurasiatherian)

test_that("conversion work as expected", {
##    skip_on_cran()
    expect_that(phy_matrix, is_a("matrix"))
    expect_that(phy_df, is_a("data.frame"))
    expect_that(phy_dnabin, is_a("DNAbin"))
    expect_that(phy_align, is_a("alignment"))
    expect_that(as.phyDat(phy_matrix), is_a("phyDat"))
    expect_that(as.phyDat(phy_df), is_a("phyDat"))
    expect_that(as.phyDat(phy_vec), is_a("phyDat"))
    expect_that(as.phyDat(phy_dnabin), is_a("phyDat"))
    expect_that(as.phyDat(phy_align), is_a("phyDat"))
    expect_is(MA_AA <- as.MultipleAlignment(chloroplast), "AAMultipleAlignment")
    expect_equal(as.phyDat(MA_AA), chloroplast)
    expect_is(MA_DNA <- as.MultipleAlignment(Laurasiatherian), "DNAMultipleAlignment")
    expect_equal(as.phyDat(MA_DNA), Laurasiatherian)
    expect_is(c2d <- codon2dna(codon_align), "phyDat")
    expect_equal(dna2codon(c2d), codon_align) 
})



test_that("subsetting and combining work as expected", {
##    skip_on_cran()
    
    expect_is(subset_1 <- subset(Laurasiatherian, select = 1:1000, site.pattern = FALSE), "phyDat")
    expect_is(subset_2 <- subset(Laurasiatherian, select = 1001:3179, site.pattern = FALSE), "phyDat")
    expect_is(lauraCbind1 <- cbind(subset_1, subset_2), "phyDat")
    expect_equal(baseFreq(lauraCbind1), baseFreq(Laurasiatherian))
    
    expect_is(subset_3 <- subset(Laurasiatherian, select = 1:100), "phyDat")
    expect_is(subset_4 <- subset(Laurasiatherian, select = 101:1605), "phyDat")
    expect_is(lauraCbind2 <- cbind(subset_3, subset_4), "phyDat")
    expect_equal(baseFreq(lauraCbind2), baseFreq(Laurasiatherian))
})


test_that("read and write work as expected", {
    skip_on_cran()
    
    write.phyDat(Laurasiatherian, "tmp1.txt")
    expect_is(laura <- read.phyDat("tmp1.txt"), "phyDat")
    expect_equal(laura, Laurasiatherian)
    unlink("tmp1.txt")
    write.phyDat(chloroplast, "tmp2.txt")
    expect_is(chloro <- read.phyDat("tmp2.txt", type="AA"), "phyDat")
    expect_equal(chloro, chloroplast)
    unlink("tmp2.txt")
    write.phyDat(chloroplast, "tmp.fas", format="fasta")
    expect_is(chloro_fas <- read.phyDat("tmp.fas", type="AA", format = "fasta"), "phyDat")
    expect_equal(chloro_fas, chloroplast)
    unlink("tmp.fas")
})

