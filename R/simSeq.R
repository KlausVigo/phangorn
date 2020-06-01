#' Simulate sequences.
#'
#' Simulate sequences from a given evolutionary tree.
#'
#' \code{simSeq} is a generic function to simulate sequence alignments
#' along a phylogeny. It is quite flexible and can generate DNA, RNA,
#' amino acids, codon, morphological or binary sequences.
#' simSeq can take as input a phylogenetic tree of class \code{phylo},
#' or a \code{pml} object; it will return an object of class \code{phyDat}.
#' There is also a more low level
#' version, which lacks rate variation, but one can combine different
#' alignments with their own rates (see example). The rate parameter acts like
#' a scaler for the edge lengths.
#'
#' For codon models \code{type="CODON"}, two additional arguments \code{dnds}
#' for the dN/dS ratio and \code{tstv} for the transition transversion ratio
#' can be supplied.
#'
#' \strong{Defaults:}
#'
#' If \code{x} is a tree of class \code{phylo}, then sequences will be generated
#' with the default Jukes-Cantor DNA model (\code{"JC"}).
#'
#' If \code{bf} is not specified, then all states will be treated as equally
#' probable.
#'
#' If \code{Q} is not specified, then a uniform rate matrix will be employed.
#'
#'
#' @param x a phylogenetic tree \code{tree}, i.e. an object of class
#' \code{phylo} or and object of class \code{pml}.
#' @param l The length of the sequence to simulate.
#' @param Q The rate matrix.
#' @param bf Base frequencies.
#' @param rootseq A vector of length \code{l} containing the root sequence.
#' If not provided, the root sequence is randomly generated.
#' @param type Type of sequences ("DNA", "AA", "CODON" or "USER").
#' @param model Amino acid model of evolution to employ, for example "WAG",
#' "JTT", "Dayhoff" or "LG". For a full list of supported models, type
#' \code{phangorn:::.aamodels}. Ignored if type is not equal to "AA".
#' @param levels A character vector of the different character tokens.
#' Ignored unless type = "USER".
#' @param rate A numerical value greater than zero giving the mutation rate
#' or scaler for edge lengths.
#' @param ancestral Logical specifying whether to return ancestral sequences.
#' @param code	The ncbi genetic code number for translation (see details). By
#' default the standard genetic code is used.
#' @param \dots Further arguments passed to or from other methods.

#' @return \code{simSeq} returns an object of class phyDat.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{phyDat}}, \code{\link{pml}}, \code{\link{SOWH.test}}
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(Laurasiatherian)
#' tree <- nj(dist.ml(Laurasiatherian))
#' fit <- pml(tree, Laurasiatherian, k=4)
#' fit <- optim.pml(fit, optNni=TRUE, model="GTR", optGamma=TRUE)
#' data <- simSeq(fit)
#' }
#'
#' tree <- rtree(5)
#' plot(tree)
#' nodelabels()
#'
#' # Example for simple DNA alignment
#' data <- simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6)
#' as.character(data)
#'
#' # Example to simulate discrete Gamma rate variation
#' rates <- discrete.gamma(1,4)
#' data1 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[1])
#' data2 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[2])
#' data3 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[3])
#' data4 <- simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[4])
#' data <- c(data1,data2, data3, data4)
#'
#' write.phyDat(data, file="temp.dat", format="sequential", nbcol = -1,
#'   colsep = "")
#' unlink("temp.dat")
#'
#' @rdname simSeq
#' @export simSeq
simSeq <- function(x, ...)
  UseMethod("simSeq")


#' @rdname simSeq
#' @method simSeq phylo
#' @export
simSeq.phylo <- function(x, l = 1000, Q = NULL, bf = NULL, rootseq = NULL,
                         type = "DNA", model = NULL, levels = NULL, rate = 1,
                         ancestral = FALSE, code=1, ...) {
  if (!is.null(model)) {
    model <- match.arg(model, .aamodels)
    getModelAA(model, bf = is.null(bf), Q = is.null(Q))
    type <- "AA"
  }

  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(extras)) {
    tmp <- c("dnds", "tstv")
    names(extras) <- tmp[pmatch(names(extras), tmp)]
    existing <- match(tmp, names(extras))
  }

  pt <- match.arg(type, c("DNA", "AA", "USER", "CODON"))
  if (pt == "DNA"){
    levels <- c("a", "c", "g", "t")
    if (!is.null(extras) ) {
      if (!is.na(existing[2]) & is.null(Q))
        tstv <- eval(extras[[existing[2]]], parent.frame())
        Q <- c(1, tstv, 1, 1, tstv, 1)
    }
  }
  if (pt == "AA")
    levels <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                "M", "F", "P", "S", "T", "W", "Y", "V")
#    c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i",
#      "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
  if (pt == "CODON") {
    .syn <- synonymous_subs(code=code)
    .sub <- tstv_subs(code=code)
    tmp <- .CODON[, as.character(code)]
    levels <- rownames(.CODON)[tmp != "*"]
    dnds <- tstv <- 1
    if (!is.null(extras)) {
      if (!is.na(existing[1]))
        dnds <- eval(extras[[existing[1]]], parent.frame())
      if (!is.na(existing[2]))
        tstv <- eval(extras[[existing[2]]], parent.frame())
    }
  }
  if (pt == "USER")
    if (is.null(levels)) stop("levels have to be supplied if type is USER")

  lbf <- length(levels)
  if (is.null(bf)) bf <- rep(1 / lbf, lbf)
  if (is.null(Q)) {
    if (type == "CODON") Q <- CodonQ(subs = .sub, syn = .syn, tstv = tstv,
                                     dnds = dnds)
    else Q <- rep(1, lbf * (lbf - 1) / 2)
  }
  if (is.matrix(Q)) Q <- Q[lower.tri(Q)]
  eig <- edQt(Q, bf)

  m <- length(levels)

  if (is.null(rootseq)) rootseq <- sample(levels, l, replace = TRUE, prob = bf)
  x <- reorder(x)
  edge <- x$edge
  nNodes <- max(edge)
  res <- matrix(NA, nNodes, l)
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  root <- as.integer(parent[!match(parent, child, 0)][1])
  res[root, ] <- rootseq
  tl <- x$edge.length
  for (i in seq_along(tl)) {
    from <- parent[i]
    to <- child[i]
    P <- getP(tl[i], eig, rate)[[1]]
    # avoid numerical problems for larger P and small t
    if (any(P < 0)) P[P < 0] <- 0
    for (j in 1:m) {
      ind <- res[from, ] == levels[j]
      res[to, ind] <- sample(levels, sum(ind), replace = TRUE, prob = P[, j])
    }
  }
  k <- length(x$tip.label)
  label <- c(x$tip.label, as.character( (k + 1):nNodes))
  rownames(res) <- label
  if (!ancestral) res <- res[x$tip.label, , drop = FALSE]
  if (pt == "DNA") return(phyDat.DNA(res, return.index = TRUE))
  if (pt == "AA") return(phyDat.AA(res, return.index = TRUE))
  if (pt == "USER") return(phyDat.default(res, levels = levels,
                                          return.index = TRUE))
  if (pt == "CODON") {
    res <- t(apply(res, 1, function(x) unlist(strsplit(x, ""))))
    return(phyDat.codon(res))
  }
}


#' @rdname simSeq
#' @method simSeq pml
#' @export
simSeq.pml <- function(x, ancestral = FALSE, ...) {
  g <- x$g
  w <- x$w
  if (x$inv > 0) {
    w <- c(x$inv, w)
    g <- c(0.0, g)
  }
  n <- length(w)
  res <- vector("list", n)
  y <- sample(n, sum(x$weight), replace = TRUE, prob = w)
  levels <- attr(x$data, "levels")
  type <- attr(x$data, "type")
  for (i in 1:n) {
    l <- sum(y == i)
    res[[i]] <- simSeq(x$tree, l, Q = x$Q, bf = x$bf, type = type,
                       levels = levels, rate = g[i], ancestral = ancestral)
  }
  x <- call("c.phyDat", quote(res[[1]]))
  if (n > 1) x <- parse(text = paste("c(", "res[[1]]",
                              paste0(",res[[", 2:n, "]]", collapse = ""), ")"))
  eval(x)
}
