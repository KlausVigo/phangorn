Sys.setenv(LANGUAGE = "en") # Force locale
library(tinytest)
using("tinysnapshot")

set.seed(123)
tree <- read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
tree2 <- read.tree(text = "((t1:1.1,t2:1.1):1,(t3:1,t4:1):1.1);")
dat <- matrix(c("a", "a",
                "a", "t",
                "t", "a",
                "t", "t"), byrow = TRUE, nrow = 4L,
              dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
dna <- phyDat(dat)
fit <- pml(tree, dna)

pml_plot <- function() plot(fit)
expect_snapshot_plot(pml_plot, "pml_plot")

trees <- c(tree, tree2) |> .compressTipLabel()

densi_plot <- function() densiTree(trees, type="phylogram", width=2,
          jitter=list(amount=0.1, random=FALSE), alpha=1)
expect_snapshot_plot(densi_plot, "densi_plot")

data("Laurasiatherian")
tmp <- subset(Laurasiatherian, 1:15, compress = TRUE)
fit <- pml_bb(tmp, "JC+G(4)", rearrangement = "NNI",
              control=pml.control(trace=0))

rates_plot <- function() plotRates(fit)
expect_snapshot_plot(rates_plot, "rates_plot")

data(woodmouse)
set.seed(123)
fit20 <- pml_bb(woodmouse, "JC", control=pml.control(trace=0),
                 ratchet.par = ratchet.control(maxit=20, bs=10))

terraces_plot <- function() terraces(fit20, pkg="plot3D")
expect_snapshot_plot(terraces_plot, "terraces_plot")
