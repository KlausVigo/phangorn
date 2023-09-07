## zzz.R

#' @import methods
#' @import Rcpp
#' @import parallel
#' @import ape
#' @importFrom stats AIC BIC logLik reorder update optim optimize constrOptim
#' @importFrom stats cophenetic hclust as.dist pchisq reshape qgamma pgamma
#' @importFrom stats na.omit model.matrix aggregate lm.fit xtabs quantile sd
#' @importFrom stats runif qbeta
#' @importFrom graphics plot plot.default plot.new plot.window text par abline
#' @importFrom graphics strwidth axis title segments points image matplot legend
#' @importFrom graphics hist identify locator barplot
#' @importFrom utils read.table download.file stack
#' @importFrom utils installed.packages write.table combn packageDescription
#' @importFrom grDevices rgb adjustcolor col2rgb
#' @useDynLib phangorn, .registration = TRUE

.packageName <- "phangorn"

.aamodels <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt",
               "MtZoa", "mtREV24", "VT", "RtREV", "HIVw", "HIVb", "FLU",
               "Blosum62", "Dayhoff_DCMut", "JTT_DCMut")

.dnamodels <- c("JC", "F81", "K80", "HKY", "TrNe", "TrN",
  "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u",
  "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe",
  "TVM", "SYM", "GTR")

.usermodels <- c("ER", "SYM", "FREQ", "GTR", "ORDERED")


# environment variables
# .CodonAlphabet <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act",
#                    "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att",
#                    "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct",
#                    "cga", "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt",
#                    "gaa", "gac", "gag", "gat", "gca", "gcc", "gcg", "gct",
#                    "gga", "ggc", "ggg", "ggt", "gta", "gtc", "gtg", "gtt",
#                    "tac", "tat", "tca", "tcc", "tcg", "tct", "tgc", "tgg",
#                    "tgt", "tta", "ttc", "ttg", "ttt")

.nucleotideAlphabet <- c("a", "c", "g", "t")

Ape_NT <- list(properties = list(
  a="a", g="g", c="c", t="t", n="n", "-"="-"),
  color=c("red", "yellow", "green", "blue", "grey", "black"))


RY_NT <- list(properties = list(
  Purine = c("a", "g", "r"),
  Pyrimidine = c("c", "t", "y"),
  "n" = "n",
  "-" = "-"),
  color=c("#FF00FF", "#00FFFF", "grey", "black"))

Ape_AA <- list(properties = list(
  Hydrophobic = c("V", "I", "L", "F", "W", "Y", "M"),
  Small = c("P", "G", "A", "C"),
  Hydrophilic = c("S", "T", "H", "N", "Q", "D", "E", "K", "R")),
  color=c("red", "yellow", "blue"))

# Properties + Conservation (Clustal X)
Clustal <- list(properties = list(
  Hydrophobic = c("A", "I", "L", "M", "F", "W", "V"),
  Positive = c("K", "R"),
  Negative = c("E", "D"),
  Polar = c("N", "Q", "S", "T"),
  Glycines = "G",
  Prolines = "P",
  Aromatic = c("H", "Y"),
  Cysteine = "C"),
  color= c("#80a0f0", "#f01505", "#c048c0", "#15c015", "#f09048", "#c0c000",
           "#15a4a4", "#f08080")
)


Polarity <- list(properties = list(
  "Non polar" = c("G", "A", "V", "L", "I", "F", "W", "M", "P"),
  "Polar, uncharged" = c("S", "T", "C", "Y", "N", "Q"),
  "Polar, acidic" = c("D", "E"),
  "Polar, basic" = c("K", "R", "H")),
  color = c("yellow", "green", "red", "blue"))

# Physicochemical Properties
Zappo_AA <- list(properties = list(
  "Aliphatic/Hydrophobic" = c("I", "L", "V", "A", "M"),
  Aromatic = c("F", "W", "Y"),
  Positive = c("K", "R", "H"),
  Negative = c("E", "D"),
  Hydrophilic = c("S", "T", "N", "Q"),
  "Conformationally special" = c("P", "G"),
  Cysteine = "C"),
  color= c("#ff7979", "#f89f56", "#0070c0", "#c00000", "#08c81a", "#cc00cc",
           "#ffff00")
)


Transmembrane_tendency <- list(properties = list(
  Lys = "K",
  Asp = "D",
  Glu = "E",
  Arg = "R",
  Gln = "Q",
  Asn = "N",
  Pro = "P",
  His = "H",
  Ser = "S",
  Thr = "T",
  Cys = "C",
  Gly = "G",
  Ala = "A",
  Tyr = "Y",
  Met = "M",
  Val = "V",
  Trp = "W",
  Leu = "L",
  Ile = "I",
  Phe = "F"
), color=c("#0000FF", "#0D00F1", "#1A00E4", "#2800D6", "#3500C9", "#4300BB",
           "#5000AE", "#5D00A1", "#6B0093", "#780086", "#860078", "#93006B",
           "#A1005D", "#AE0050", "#BB0043", "#C90035", "#D60028", "#E4001A",
           "#F1000D", "#FF0000"))


# if rate g[i] is smaller than .gEps invariant site is increased by w[i]
.gEps <- 1e-12

.PlotNetworxEnv <- new.env()


loadModule("Fitch_mod", TRUE)

# .onLoad  <- function(libname, pkgname) {
#    library.dynam("phangorn", pkgname, libname)
#}

