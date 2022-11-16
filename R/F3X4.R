# compute frequencies for codon positions
# TODO: using baseFreq codon should save memory instead of codon2dna
bf_by_codon <- function(x) {
  y <- codon2dna(x)
  ny <- sum(attr(y, "weight"))
  M <- matrix(0, nrow = 4, ncol = 3)
  M[, 1] <- subset(y, select = seq(1, ny, by = 3), site.patter = FALSE) |>
    baseFreq()
  M[, 2] <- subset(y, select = seq(2, ny, by = 3), site.patter = FALSE) |>
    baseFreq()
  M[, 3] <- subset(y, select = seq(3, ny, by = 3), site.patter = FALSE) |>
    baseFreq()
  M
}

# return frequencies for all 61 states
F3x4_freq <- function(M, CodonAlphabet,
                      nucleotideAlphabet = .nucleotideAlphabet) {
  pos <- CodonAlphabet |> strsplit("") |> unlist() |>
    match(nucleotideAlphabet) |> matrix(ncol = 3, byrow = TRUE)
  codon_frequencies <- M[pos[, 1], 1] * M[pos[, 2], 2] * M[pos[, 3], 3]
  codon_frequencies / sum(codon_frequencies)
}


F3x4 <- function(x) {
  BF <- bf_by_codon(x)
  codon_abc <- attr(x, "levels")
  F3x4_freq(BF, CodonAlphabet = codon_abc)
}


F1x4 <- function(x, bf=NULL) {
  if(is.null(bf))bf <- codon2dna(x) |> baseFreq()
  BF <- matrix(bf, 4, 3)
  codon_abc <- attr(x, "levels")
  F3x4_freq(BF, CodonAlphabet = codon_abc)
}
