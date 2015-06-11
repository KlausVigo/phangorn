context("conversion_and_subsetting")

data(Laurasiatherian)
phy_matrix <- as.character(Laurasiatherian)
phy_df <- as.data.frame(Laurasiatherian)
phy_dnabin <- as.DNAbin(Laurasiatherian)
phy_align <- phyDat2alignment(Laurasiatherian)

test_that("conversion work as expected", {
    skip_on_cran()
    expect_that(phy_matrix, is_a("matrix"))
    expect_that(phy_df, is_a("data.frame"))
    expect_that(phy_dnabin, is_a("DNAbin"))
    expect_that(phy_align, is_a("alignment"))
    expect_that(as.phyDat(phy_matrix), is_a("phyDat"))
    expect_that(as.phyDat(phy_df), is_a("phyDat"))
    expect_that(as.phyDat(phy_dnabin), is_a("phyDat"))
    expect_that(as.phyDat(phy_align), is_a("phyDat"))
})