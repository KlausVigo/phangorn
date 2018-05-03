## zzz.R 

.packageName <- "phangorn"

.aamodels <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU","Blosum62","Dayhoff_DCMut","JTT_DCMut")

.dnamodels <- c("JC", "F81", "K80", "HKY", "TrNe", "TrN", 
  "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", 
  "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", 
  "TVM", "SYM", "GTR")


# if rate g[i] is smaller than .gEps invariant site is increased by w[i]
.gEps <- 1e-12


.PlotNetworxEnv <- new.env()

# .onLoad  <- function(libname, pkgname) {
#    library.dynam("phangorn", pkgname, libname)
#}
