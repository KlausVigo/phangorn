# branch and bound
getOrder <- function(x) {
  label <- names(x)
  dm <- as.matrix(dist.hamming(x, FALSE))
  ind <- as.vector(which(dm == max(dm), arr.ind = TRUE)[1, ])
  nTips <- as.integer(length(label))
  added <- ind
  remaining <- c(1:nTips)[-ind]

  tree <- structure(list(edge = structure(c(rep(nTips + 1L, 3), c(ind, 0L)),
        .Dim = c(3L, 2L)), tip.label = label, Nnode = 1L), .Names = c("edge",
        "tip.label", "Nnode"), class = "phylo", order = "postorder")

  l <- length(remaining)
  res <- numeric(l)

  nr <- attr(x, "nr")
  storage.mode(nr) <- "integer"

  weight <- attr(x, "weight")
  storage.mode(weight) <- "double"

  f <- init_fitch(x, FALSE, FALSE, m=4L)

  edge <- tree$edge
  for (i in seq_along(remaining)) {
    edge[3, 2] <- remaining[i]
    res[i] <- f$pscore(edge)
  }
  tmp <- which.max(res)
  added <- c(added, remaining[tmp])
  remaining <- remaining[-tmp]
  tree$edge[, 2] <- added

  while (length(remaining) > 0) {
    edge <- tree$edge[, 2] + 2 * nTips
    f$prep_spr(tree$edge)
    l <- length(remaining)
    res <- numeric(l)
    nt <- numeric(l)
    for (j in 1:l) {
      score <- f$pscore_vec(edge, remaining[j])
      res[j] <- min(score)
      nt[j] <- which.min(score)
    }
    tmp <- which.max(res)
    added <- c(added, remaining[tmp])
    tree <- addOne(tree, remaining[tmp], nt[tmp])
    remaining <- remaining[-tmp]
  }
  added <- c(added, remaining)
  added
}


seq_stats <- function(x){
  nr <- attr(x, "nr")
  contrast <- attr(x, "contrast")
  a <- seq_len(nr)
  STATE <- POS <- matrix(0L, nrow(contrast), nr)
  for(i in seq_along(x)){
    IND <- cbind(x[[i]], a)
    STATE[IND] <- STATE[IND] + 1L
    POS[IND] <- i
  }
  list(state=STATE, position=POS)
}

# Incompatibility lower Bound
ilb <- function(x, LB) {
  nr <- attr(x, "nr")
  contrast <- attr(x, "contrast")
  rownames(contrast) <- attr(x, "allLevels")
  colnames(contrast) <- attr(x, "levels")
  weight0 <- attr(x, "weight")
  attr(x, "weight") <- rep(1, nr)
  attr(x, "index") <- NULL

  y <- as.character(x)
  singles <- attr(x, "levels")
  fun2 <- function(x, singles) all(x %in% singles)
  fun1 <- function(x) cumsum(!duplicated(x)) - 1L
  fun3 <- function(x) sum(!duplicated(x)) - 1L

  tmp <- apply(y, 2, fun2, singles)
  ind <- which(tmp)
  if (length(ind) < 2) return(0)

  y <- y[, ind, drop = FALSE]
  weight0 <- weight0[ind]
  single_dis <- LB[, ind]
  nTips <- nrow(y)
  l <- length(weight0)
  res <- numeric(nTips)
  for (i in 1:(l - 1)) {
    for (j in (i + 1):l) {
      if ((weight0[i] > 0) & (weight0[j] > 0)) {
        z <- paste(y[, i], y[, j], sep = "_")
        dis2 <- single_dis[, i] + single_dis[, j]
        dis <- fun1(z)
        if (dis[nTips] > dis2[nTips]) {
          dis <- pmax(dis, dis2) - dis2
          if (sum(dis[4:nTips]) > 0) {
            wmin <- min(weight0[i], weight0[j])
            weight0[i] <- weight0[i] - wmin
            weight0[j] <- weight0[j] - wmin
            res <- res + dis * wmin
          }
        }
      }
      if(weight0[i] < 1e-6) break()
    }
  }
  res
}


#' Branch and bound for finding all most parsimonious trees
#'
#' \code{bab} finds all most parsimonious trees.
#'
#' This implementation is very slow and depending on the data may take very
#' long time. In the worst case all \eqn{(2n-5)!! = 1 \times 3 \times 5
#' \times \ldots \times (2n-5)}{(2n-5)!! = 1 * 3 * 5 * ... * (2n-5)} possible
#' trees have to be examined, where n is the number of species / tips. For ten
#' species there are already 2027025 tip-labelled unrooted trees. It only uses
#' some basic strategies to find a lower and upper bounds similar to penny from
#' phylip. \code{bab} uses a very basic heuristic approach of MinMax Squeeze
#' (Holland et al. 2005) to improve the lower bound.  On the positive side
#' \code{bab} is not like many other implementations restricted to binary or
#' nucleotide data.
#'
#' @aliases bab BranchAndBound
#' @param data an object of class phyDat.
#' @param tree a phylogenetic tree an object of class phylo, otherwise a
#' pratchet search is performed.
#' @param trace defines how much information is printed during optimization.
## @param ILBound compute incompatibility lower bound (default TRUE) of
## Holland (2005).
#' @param \dots Further arguments passed to or from other methods
#' @return \code{bab} returns all most parsimonious trees in an object of class
#' \code{multiPhylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} based on work on Liam
#' Revell
#' @seealso \code{\link{pratchet}}, \code{\link{dfactorial}}
#' @references Hendy, M.D. and Penny D. (1982) Branch and bound algorithms to
#' determine minimal evolutionary trees.  \emph{Math. Biosc.} \bold{59},
#' 277-290
#'
#' Holland, B.R., Huber, K.T. Penny, D. and Moulton, V. (2005) The MinMax
#' Squeeze: Guaranteeing a Minimal Tree for Population Data, \emph{Molecular
#' Biology and Evolution}, \bold{22}, 235--242
#'
#' White, W.T. and Holland, B.R. (2011) Faster exact maximum parsimony search
#' with XMP. \emph{Bioinformatics}, \bold{27(10)},1359--1367
#' @keywords cluster
#' @examples
#'
#' data(yeast)
#' dfactorial(11)
#' # choose only the first two genes
#' gene12 <- yeast[, 1:3158]
#' trees <- bab(gene12)
#'
#' @export
bab <- function(data, tree = NULL, trace = 0, ...) {
  if (hasArg(ILBound))
    ILBound <- list(...)$ILBound
  else ILBound <- FALSE
  if(inherits(data, "DNAbin") | inherits(data, "AAbin")) data <- as.phyDat(data)
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  compress <- TRUE
  recursive <- TRUE
  nTips <- length(data)
  if (nTips < 4) return(stree(nTips, tip.label = names(data)))

  data <- removeParsimonyUninfomativeSites(data, recursive=recursive)
  star_tree <- ifelse(attr(data, "nr") == 0, TRUE, FALSE)
  add_taxa <- ifelse(is.null(attr(data, "duplicated")), FALSE, TRUE)
  p0 <- attr(data, "p0")

  nTips <- length(data)
  if (nTips < 4L  || star_tree) {
    nam <- names(data)
    if (star_tree) tree <- stree(length(nam), tip.label = nam)
    else tree <- stree(nTips, tip.label = names(data))
    if(add_taxa) tree <- addTaxa(tree, attr(data, "duplicated"))
    tree <- unroot(tree)
    return(tree)
  }

  # compress sequences (all transitions count equal)
  data <- compressSites(data)
  o <- order(attr(data, "weight"), decreasing = TRUE)
  data <- subset(data, select = o, site.pattern=TRUE)

  nr <- as.integer(attr(data, "nr"))
  inord <- getOrder(data)
  data <- data[inord,]
  tree <- pratchet(data, start = tree, trace = trace - 1, maxit=10,
                   all=FALSE, ...)
  p_vec <- fitch(tree, data, "site")

  nTips <- m <- length(data)
  TMP <- matrix(0, m, nr) # UB <-
  for (i in 2:m) {
    TMP[i, ] <- lowerBound(data[1:i,])
  }

  weight <- as.double(attr(data, "weight"))

  m <- nr * (2L * nTips - 2L)
  # Single column discrepancy
  mmsAmb <- TMP %*% weight
  #  mmsAmb <- mmsAmb[nTips] - mmsAmb
  mms0 <- mms1 <- 0
  if (ILBound) mms1 <- ilb(data, TMP)
  mms0 <- mms1 + mmsAmb
  mms0 <- mms0[nTips] - mms0
  mms0 <- c(mms0, 0)

  f <- init_fitch(data, m=4L)

  if (trace > 1) print(paste("lower bound:", p0 + mms0[1]))
  bound <- f$pscore(tree$edge)
  if (trace > 1) print(paste("upper bound:", bound + p0))

  startTree <- structure(list(edge = structure(c(rep(nTips + 1L, 3),
        as.integer(1:3)), .Dim = c(3L, 2L)), tip.label = names(data),
        Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"), class = "phylo",
        order = "postorder")

  trees <- vector("list", nTips)
  trees[[3]] <- list(startTree$edge)
  for (i in 4:nTips) trees[[i]] <- vector("list", (2L * i) - 5L) # new

  # index M[i] is neues node fuer edge i+1
  # index L[i] is length(node) tree mit i+1
  L <- as.integer(2L * (1L:nTips) - 3L)
  M <- as.integer(1L:nTips + nTips - 1L)

  PSC <- matrix(c(3, 1, 0), 1, 3)
  PSC[1, 3] <- f$pscore(startTree$edge)

  k <- 4L
  Nnode <- 1L
  npsc <- 1

  visited <- numeric(nTips)

  result <- list()
  while (npsc > 0) {
    a <- PSC[npsc, 1]
    b <- PSC[npsc, 2]
    blub <- PSC[npsc, 3]
    PSC <- PSC[-npsc, , drop = FALSE]
    npsc <- npsc - 1L
    tmpTree <- trees[[a]][[b]]
    edge <- tmpTree[, 2] + 2 * nTips

    f$prep_spr(tmpTree)
    score <- f$pscore_vec(edge, as.integer(a + 1L))
    score <- score + blub + mms0[a + 1L]
    ms <- min(score)
    if (ms < bound + .1) {
      if ((a + 1L) < nTips) {
        ind <- (1:L[a])[score <= bound]
        trees[[a + 1]][seq_along(ind)] <- .Call('AddOnes', tmpTree,
                        as.integer(a + 1L), as.integer(ind), as.integer(L[a]),
                        as.integer(M[a]))
        l <- length(ind)
        # os <- order(score[ind], decreasing=TRUE)
        os <- seq_len(l)
        # in C++ pushback
        PSC <- rbind(PSC, cbind(rep(a + 1, l), os, score[ind] - mms0[a + 1L]))
        npsc <- npsc + l
        visited[a + 1] <- visited[a + 1] + l
        #  PSC = rbind(PSC, cbind(rep(a+1, l), os, score[ind][os] ))
      }
      else {
        ind <- which(score == ms)
        tmp <- vector("list", length(ind))
        tmp[seq_along(ind)] <- .Call('AddOnes', tmpTree,
                                     as.integer(a + 1L), as.integer(ind),
                                     as.integer(L[a]), as.integer(M[a]))
        if (ms < bound) {
          bound <- ms
          if (trace) cat("upper bound:", bound + p0, "\n")
          result <- tmp
          PSC <- PSC[PSC[, 3] < (bound + 1e-8), ]
          npsc <- nrow(PSC)
        }
        else result <- c(result, tmp)
      }
    }
  }
  for (i in seq_along(result)) {
    result[[i]] <- structure(list(edge = result[[i]], Nnode = nTips - 2L),
                             .Names = c("edge", "Nnode"), class = "phylo",
                             order = "postorder")
  }
  attr(result, "TipLabel") <- tree$tip.label
  class(result) <- "multiPhylo"
  if(add_taxa) result <- addTaxa(result, attr(data, "duplicated"))
  attr(result, "visited") <- visited
  return(result)
}
