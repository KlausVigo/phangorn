Sys.setenv(LANGUAGE = "en") # Force locale
library(tinytest)
using("tinysnapshot")


tree <- read.tree(text = "(((t1:1,t2:1):1,t3:2):1,(t4:1.5, t5:1.5):1.5);")
dat <- matrix(c("a", "c",
                "a", "g",
                "t", "c",
                "t", "g",
                "t", "g"), byrow = TRUE, nrow = 5L,
              dimnames = list(c("t1", "t2", "t3", "t4", "t5"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)
test_ml <- anc_pml(fit)



pie_ancestral <- function() {
  plotAnc(test_ml, pos="bottomleft")
  add_mutations(test_ml, adj = c(0.5, -0.3))
}
expect_snapshot_plot(pie_ancestral, "pie_ancestral")



old_par <- par(mar=c(5,4,4,8))
heat <- function() anc_heatmap(test_ml)
expect_snapshot_plot(heat, "heat")
par(old_par)

#seq_logo <- function()plotSeqLogo(test_ml)
#expect_snapshot_plot(seq_logo, "seq_logo")
