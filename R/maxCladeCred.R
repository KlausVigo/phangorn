#' Maximum clade credibility tree
#'
#' \code{maxCladeCred} computes the maximum clade credibility tree from a
#' sample of trees.
#' So far just the best tree is returned. No annotations or transformations of
#' edge length are performed and the edge length are taken from the tree.
#'
#' If a list of partition is provided then the clade credibility is computed
#' for the trees in x.
#'
#' \code{allCompat} returns a 50\% majority rule consensus tree with added
#' compatible splits similar to the option allcompat in MrBayes. This tree has
#' no edge length.
#'
#' \code{\link{add_edge_length}} can be used to add edge lengths computed from
#' the sample of trees.
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
#' @seealso \code{\link[ape]{consensus}}, \code{\link{consensusNet}},
#' \code{\link[ape]{prop.part}}, \code{\link{bootstrap.pml}},
#' \code{\link{plotBS}}, \code{\link{transferBootstrap}},
#' \code{\link{add_edge_length}}, \code{\link{add_boxplot}}
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
#'
#' par(mfrow = c(2,1))
#' plot(max_clade_cred, main="Edge length from tree")
#' add_boxplot(max_clade_cred, bs)
#' max_clade_cred_2 <- add_edge_length(max_clade_cred, bs)
#' plot(max_clade_cred_2, main="Edge length computed from sample")
#' add_boxplot(max_clade_cred_2, bs)
#'
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
    tr <- addConfidences(tr, pp, rooted=rooted)
    attr(tr, "clade.credibility") <- res[k]
    return(tr)
  }
  res
}


#' @rdname maxCladeCred
#' @export
mcc <- maxCladeCred


#' @rdname maxCladeCred
#' @export
allCompat <- function(x, rooted=FALSE) {
  x <- .compressTipLabel(x)
  if(!rooted) x <- unroot(x, collapse.singles=TRUE)
  l <- length(x)
  nt <- Ntip(x[[1]])
  pp <- prop.part(x)
  if(!rooted) pp <- postprocess.prop.part(pp, method = "SHORTwise")
  spl <- as.splits(pp)
  if(rooted) spl <- SHORTwise(spl)
  ord <- order(attr(spl, "weights"), decreasing = TRUE)
  spl <- spl[ord]
  if(rooted){
    ind <- duplicated(spl) | (lengths(spl)==1)
    root_candidates <- spl[ind]
    spl <- spl[!ind]
  }
  ll <- lengths(spl)
  spl <- spl[ll>0]
  w <- attr(spl, "weights")
  ind <- (w / l) > 0.5 & (lengths(spl) > 1)
  if(!any(ind)) ind[which.min(w)] <- TRUE
  res <- spl[ind]
  spl <- spl[!ind]
  for(i in seq_along(spl)){
    if(all(compatible(res, spl[i]) == 0)) res <- c(res, spl[i])
  }
  tree <- as.phylo(res, FALSE)
  if(rooted && length(root_candidates)>0){
    labels <- attr(x, "TipLabel")[root_candidates[[1]]]
    tree <- root(tree, outgroup = labels, resolve.root = TRUE)
  }
  tree$edge.length <- NULL
  tree <- addConfidences(tree, pp, rooted=rooted)
  tree
}
