context("conversion_and_subsetting")
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

test_that("conversion work as expected", {
##    skip_on_cran()
    expect_is(phy_matrix, "matrix")
    expect_is(phy_df, "data.frame")
    expect_s3_class(phy_dnabin, "DNAbin")
    expect_s3_class(phy_aabin, "AAbin")
    expect_s3_class(phy_align, "alignment")
    expect_s3_class(as.phyDat(phy_matrix), "phyDat")
    expect_s3_class(as.phyDat(phy_df), "phyDat")
    expect_s3_class(as.phyDat(phy_dnabin), "phyDat")
    expect_equal(as.phyDat(phy_aabin), chloroplast)
    expect_s3_class(phyDat(phy_vec_dna), "phyDat")
    expect_s3_class(phyDat(phy_vec_user, type="USER", levels = c("0","1")), "phyDat")
    expect_s3_class(as.phyDat(phy_dnabin), "phyDat")
    expect_s3_class(as.phyDat(phy_align), "phyDat")
    expect_s3_class(c2d <- codon2dna(codon_align), "phyDat")
    expect_equal(dna2codon(c2d), codon_align)
})


test_that("conversion with Biostrings work as expected", {
    skip_on_cran()
    if(requireNamespace('Biostrings')){
      expect_s4_class(MA_AA <- as.MultipleAlignment(chloroplast),
                  "AAMultipleAlignment")
        expect_equal(as.phyDat(MA_AA), chloroplast)
        expect_s4_class(MA_DNA <- as.MultipleAlignment(Laurasiatherian),
                  "DNAMultipleAlignment")
        expect_equal(as.phyDat(MA_DNA), Laurasiatherian)
    }
})


test_that("subsetting and combining work as expected", {
##    skip_on_cran()

  expect_s3_class(subset_1 <- subset(Laurasiatherian, select = 1:1000,
                                 site.pattern = FALSE), "phyDat")
  expect_s3_class(subset_2 <- subset(Laurasiatherian, select = 1001:3179,
                                 site.pattern = FALSE), "phyDat")
  expect_s3_class(lauraCbind1 <- cbind(subset_1, subset_2), "phyDat")
  expect_equal(baseFreq(lauraCbind1), baseFreq(Laurasiatherian))

  expect_s3_class(subset_3 <- subset(Laurasiatherian, select = 1:100), "phyDat")
  expect_s3_class(subset_4 <- subset(Laurasiatherian, select = 101:1605), "phyDat")
  expect_equal(subset_1, Laurasiatherian[, 1:1000])
  expect_s3_class(lauraCbind2 <- cbind(subset_3, subset_4), "phyDat")
  expect_equal(baseFreq(lauraCbind2), baseFreq(Laurasiatherian))
})


test_that("read and write work as expected", {
    skip_on_cran()

    write.phyDat(Laurasiatherian, "tmp1.txt")
    expect_s3_class(laura <- read.phyDat("tmp1.txt"), "phyDat")
    expect_equal(laura, Laurasiatherian)
    unlink("tmp1.txt")
    write.phyDat(chloroplast, "tmp2.txt")
    expect_s3_class(chloro <- read.phyDat("tmp2.txt", type="AA"), "phyDat")
    expect_equal(chloro, chloroplast)
    unlink("tmp2.txt")
    write.phyDat(chloroplast, "tmp.fas", format="fasta")
    expect_s3_class(chloro_fas <- read.phyDat("tmp.fas", type="AA",
                                              format = "fasta"), "phyDat")
    expect_equal(chloro_fas, chloroplast)
    unlink("tmp.fas")
})


test_that("removing duplicated sequences works", {
    skip_on_cran()

    tmp <- as.character(Laurasiatherian)
    laura <- phyDat(rbind(phy_matrix, phy_matrix))
    names(laura) <- paste0(names(laura), rep(c(1,2), each=47))

    map1 <- map_duplicates(laura)
    map2 <- map_duplicates(Laurasiatherian)
    expect_null(map2)
    expect_is(map1, "data.frame")
})
