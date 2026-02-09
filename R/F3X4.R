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


# pml, update.pml
# bf_choice <- match.arg(j, c("equal", "empirical", "F1x4", "F3x4", "F61"))
bf_char <- function(data, bf_choice ){
  nc <- attr(data, "nc")
  type <- attr(data, "type")
  txt <-  deparse(substitute(bf_choice))
  if(bf_choice %in% c("F1x4", "F3x4", "F61") && type != "CODON")
    stop(gettextf("%s not available for this data type", txt))
  bf <- switch(bf_choice,
               equal = rep(1 / nc, nc),
               empirical = baseFreq(data),
               F61 = baseFreq(data),
               F3x4 = F3x4(data),
               F1x4 = F1x4(data))
  if(has_gap_state(data) && bf_choice=="equal"){
    bf <- baseFreq(data)
    bf[-nc] <- (1 - bf[nc]) / (nc-1)
  }
  names(bf) <- NULL
  bf
}
# freq_df <- df_freq_codon(bf_choice)

#
updateModel <- function(model, type, nc, bf = TRUE, Q = TRUE,
                        has_gap_state = FALSE, initiate = FALSE){
  if (type == "AA") {
    model <- match.arg(eval(model), .aa_3Di_models)
    tmp <- get(paste0(".", model), environment(pml))
    if (has_gap_state) {
      tmp$Q <- add_gap_Q_AA(tmp$Q)
      tmp$bf <- add_gap_bf_AA(tmp$bf)
    }
    if (Q) assign("Q", tmp$Q, envir = parent.frame())
    if (bf) assign("bf", tmp$bf, envir = parent.frame())
    return(NULL)
  } else if (type %in% c("DNA", "USER")) {
    if (type == "USER") {
      model <- match.arg(model, .usermodels)
      sc <- subsChoice_USER(model, nc)
    } else{
      model <- match.arg(model, .dnamodels)
      sc <- subsChoice(model, has_gap_state)
    }
    tmpbf <- ifelse(sc$optBf, "empirical", "equal")
    tmpQ <- rep(1, (nc * (nc - 1L)/2))
#    browser()
    if (model == "ORDERED") tmpQ <- sc$Q
    if (Q && !sc$optQ) assign("Q", tmpQ, envir = parent.frame())
    if(initiate){
      if (bf && sc$optBf) assign("bf", tmpbf, envir = parent.frame())
    } else {
      if (bf && !sc$optBf) assign("bf", tmpbf, envir = parent.frame())
    }

  }
}

