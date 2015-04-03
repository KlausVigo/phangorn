## zzz.R 

.packageName <- "phangorn"

.aamodels <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU","Blosum62","Dayhoff_DCMut","JTT_DCMut")


# if g[i] is smaller .gEps inv is increased w[i]
.gEps <- 1e-30


# .onLoad  <- function(libname, pkgname) {
#    library.dynam("phangorn", pkgname, libname)
#}
