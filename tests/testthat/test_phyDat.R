context("conversion_and_subsetting")

data(Laurasiatherian)
data(chloroplast)
phy_matrix <- as.character(Laurasiatherian)
phy_df <- as.data.frame(Laurasiatherian)
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
    expect_that(as.phyDat(phy_dnabin), is_a("phyDat"))
    expect_that(as.phyDat(phy_align), is_a("phyDat"))
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
})

