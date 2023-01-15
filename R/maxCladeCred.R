#' Maximum clade credibility tree
#'
#' \code{maxCladeCred} computes the maximum clade credibility tree from a
#' sample of trees.
#'
#' So far just the best tree is returned. No annotations or transformations of
#' edge length are performed.
#'
#' If a list of partition is provided then the clade credibility is computed
#' for the trees in x.
#'
#' \code{allCompat} returns a 50\% majority rule consensus tree with added
#' compatible splits similar to the option allcompat in MrBayes.
#'
#' @param x \code{x} is an object of class \code{multiPhylo} or \code{phylo}
#' @param tree logical indicating whether return the tree with the clade
#' credibility (default) or the clade credibility score for all trees.
#' @param rooted logical, if FALSE the tree with highest maximum bipartition
#' credibility is returned.
#' @param part a list of partitions as returned by \code{prop.part}
#' @return a tree (an object of class \code{phylo}) with the highest clade
#' credibility or a numeric vector of clade credibilities for each tree.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{consensus}}, \code{\link{consensusNet}},
#' \code{\link{prop.part}}, \code{\link{bootstrap.pml}}, \code{\link{plotBS}},
#' \code{\link{transferBootstrap}}
#' @keywords cluster
#' @importFrom fastmatch fmatch
#' @examples
#'
#'
#' data(Laurasiatherian)
#' set.seed(42)
#' bs <- bootstrap.phyDat(Laurasiatherian,
#'   FUN = function(x)upgma(dist.hamming(x)), bs=100)
#'
#' strict_consensus <- consensus(bs)
#' majority_consensus <- consensus(bs, p=.5)
#' all_compat <- allCompat(bs)
#' max_clade_cred <- maxCladeCred(bs)
#'
#' old.par <- par(no.readonly = TRUE)
#' par(mfrow = c(2,2), mar = c(1,4,1,1))
#' plot(strict_consensus, main="Strict consensus tree")
#' plot(majority_consensus, main="Majority consensus tree")
#' plot(all_compat, main="Majority consensus tree with compatible splits")
#' plot(max_clade_cred, main="Maximum clade credibility tree")
#' par(old.par)
#'
#' # compute clade credibility for trees given a prop.part object
#' pp <- prop.part(bs)
#' tree <- rNNI(bs[[1]], 20)
#' maxCladeCred(c(tree, bs[[1]]), tree=FALSE, part = pp)
#' # first value likely be -Inf
#'
#' @rdname maxCladeCred
#' @export
maxCladeCred <- function(x, tree = TRUE, part = NULL, rooted = TRUE) {
  if (inherits(x, "phylo")) x <- c(x)
  if (is.null(part)) {
    if (!rooted) {
      pp <- unroot(x) |> prop.part()
    } else {
      pp <- prop.part(x)
    }
  }
  else {
    pp <- part
  }
  pplabel <- attr(pp, "labels")
  if (!rooted){
    pp <- postprocess.prop.part(pp, method="SHORTwise")
  }
  x <- .uncompressTipLabel(x)
  class(x) <- NULL
  m <- max(attr(pp, "number"))
  nb <- log(attr(pp, "number") / m)
  l <- length(x)
  res <- numeric(l)
  for (i in 1:l) {
    tmp <- checkLabels(x[[i]], pplabel)
    if (!rooted) tmp <- unroot(tmp)
    ppi <- prop.part(tmp) # trees[[i]]
    if (!rooted) ppi <- SHORTwise(ppi)
    indi <- fmatch(ppi, pp)
    if (any(is.na(indi))) {
      res[i] <- -Inf
    } else {
      res[i] <- sum(nb[indi])
    }
  }
  if (tree) {
    k <- which.max(res)
    tr <- x[[k]]
    tr <- addConfidences(tr, pp)
    attr(tr, "clade.credibility") <- res[k]
    return(tr)
  }
  res
}


#' @rdname maxCladeCred
#' @export
mcc <- maxCladeCred


compatible_clades <- function(obj1, obj2) {
  if (!inherits(obj1, "splits"))
    stop("obj needs to be of class splits")
  labels <- attr(obj1, "labels")
  l <- length(labels)
  n <- length(obj1)
  bp1 <- as.matrix(obj1)

  rs1 <- rowSums(bp1)

  m <- length(obj2)
  bp2 <- as.matrix(obj2)
  labels2 <- attr(obj2, "labels")
  bp2 <- bp2[, match(labels2, labels), drop = FALSE]

  rs2 <- rowSums(bp2)

  bp3 <- bp2
  bp3[bp3[, 1] == 0L, ] <- 1L - bp3[bp3[, 1] == 0L, ]

  res <- matrix(0L, n, m)

  tmp1 <- tcrossprod(bp1, bp2)
  tmp2 <- as.integer(!(tmp1==rs1))
  tmp3 <- as.integer(!(tmp1==rs2))
  res[(tmp1 * tmp2 * tmp3) > 0] <- 1L
  return(res)
}


#' @rdname maxCladeCred
#' @export
allCompat <- function(x, rooted=FALSE) {
  x <- unroot(x)
  l <- length(x)
  pp <- prop.part(x)
  if(!rooted) pp <- postprocess.prop.part(pp, method = "SHORTwise")
  spl <- as.splits(pp)
  w <- attr(spl, "weights")
  ind <- (w / l) > 0.5
  res <- spl[ind]
  spl <- spl[!ind]
  w <- attr(spl, "weights")
  ord <- order(w, decreasing = TRUE)
  for(i in ord){
    if(rooted){
      if(all(compatible_clades(res, spl[i]) == 0)) res <- c(res, spl[i])
    }
    else if(all(compatible(res, spl[i]) == 0)) res <- c(res, spl[i])
  }
  tree <- as.phylo(res, FALSE)
  if(rooted){
    ll <- lengths(res)
    res <- res[ll < length(attr(res, "labels"))]
    ind <- which.max(lengths(res))
    tree <- root(tree, outgroup = res[[ind]], resolve.root = TRUE)
  }
  tree$edge.length <- NULL
  tree <- addConfidences(tree, pp)
  tree
}
