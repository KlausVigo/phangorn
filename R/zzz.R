## zzz.R

#' @useDynLib phangorn, .registration = TRUE

.packageName <- "phangorn"

.aamodels <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt",
               "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU",
               "Blosum62","Dayhoff_DCMut","JTT_DCMut")

.dnamodels <- c("JC", "F81", "K80", "HKY", "TrNe", "TrN",
  "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u",
  "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe",
  "TVM", "SYM", "GTR")

# environment variables
.CodonAlphabet <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act",
                    "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att",
                    "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct",
                    "cga", "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt",
                    "gaa", "gac", "gag", "gat", "gca", "gcc", "gcg", "gct",
                    "gga", "ggc", "ggg", "ggt", "gta", "gtc", "gtg", "gtt",
                    "tac", "tat", "tca", "tcc", "tcg", "tct", "tgc", "tgg",
                    "tgt", "tta", "ttc", "ttg", "ttt")

.nucleotideAlphabet <- c("a", "c", "g", "t")



# if rate g[i] is smaller than .gEps invariant site is increased by w[i]
.gEps <- 1e-12

.PlotNetworxEnv <- new.env()


# loadModule("Fitch_mod", TRUE)

# .onLoad  <- function(libname, pkgname) {
#    library.dynam("phangorn", pkgname, libname)
#}
