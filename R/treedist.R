#
# tree distance functions
#
coph <- function(x, path = FALSE) {
  if (is.null(attr(x, "order")) || attr(x, "order") != "postorder")
    x <- reorder(x, "postorder")
  el <- x$edge.length
  if (path) el <- rep(1.0, nrow(x$edge))
  nTips <- as.integer(length(x$tip.label))
  nNode <- as.integer(x$Nnode)
  dm <- cophenetic_cpp(x$edge, as.double(el), nTips, nNode)
  attr(dm, "Size") <- nTips
  attr(dm, "Labels") <- x$tip.label
  attr(dm, "Diag") <- FALSE
  attr(dm, "Upper") <- FALSE
  class(dm) <- "dist"
  dm
}


#' @export
cophenetic.splits <- function(x) {
  labels <- attr(x, "labels")
  X <- splits2design(x)
  weights <- attr(x, "weight")
  if(is.null(weights)) weights <- rep(1, length(x))
  dm <- as.vector(X %*% weights)
  attr(dm, "Size") <- length(labels)
  attr(dm, "Labels") <- labels
  attr(dm, "Diag") <- FALSE
  attr(dm, "Upper") <- FALSE
  class(dm) <- "dist"
  dm
}



#' Pairwise Distances from a Phylogenetic Network
#'
#' \code{cophenetic.networx} computes the pairwise distances between the pairs
#' of tips from a phylogenetic network using its branch lengths.
#'
#'
#' @aliases cophenetic.networx cophenetic.splits
#' @param x an object of class \code{networx}.
#' @return an object of class \code{dist}, names are set according to the tip
#' labels (as given by the element \code{tip.label} of the argument \code{x}).
#' @author Klaus Schliep
#' @seealso \code{\link[stats]{cophenetic}} for the generic function,
#' \code{neighborNet} to construct a network from a distance matrix
#' @keywords manip
#' @export
cophenetic.networx <- function(x) {
  spl <- x$splits
  cophenetic.splits(spl)
}


## @aliases treedist RF.dist wRF.dist KF.dist path.dist sprdist SPR.dist
#' Distances between trees
#'
#' \code{treedist} computes different tree distance methods and \code{RF.dist}
#' the Robinson-Foulds or symmetric distance. The Robinson-Foulds distance only
#' depends on the topology of the trees. If edge weights should be considered
#' \code{wRF.dist} calculates the weighted RF distance (Robinson & Foulds
#' 1981). and \code{KF.dist} calculates the branch score distance (Kuhner &
#' Felsenstein 1994).  \code{path.dist} computes the path difference metric as
#' described in Steel and Penny 1993).
#' \code{sprdist} computes the approximate SPR distance (Oliveira Martins et
#' al. 2008, de Oliveira Martins 2016).
#'
#' @details The Robinson-Foulds distance between two trees \eqn{T_1} and \eqn{T_2} with
#' \eqn{n} tips is defined as (following the notation Steel and Penny 1993):
#' \deqn{d(T_1, T_2) = i(T_1) + i(T_2) - 2v_s(T_1, T_2)} where \eqn{i(T_1)}
#' denotes the number of internal edges and \eqn{v_s(T_1, T_2)} denotes the
#' number of internal splits shared by the two trees. The normalized
#' Robinson-Foulds distance is derived by dividing \eqn{d(T_1, T_2)} by the
#' maximal possible distance \eqn{i(T_1) + i(T_2)}. If both trees are unrooted
#' and binary this value is \eqn{2n-6}.
#'
#' Functions like \code{RF.dist} returns the Robinson-Foulds distance (Robinson
#' and Foulds 1981) between either 2 trees or computes a matrix of all pairwise
#' distances if a \code{multiPhylo} object is given.
#'
#' For large number of trees the distance functions can use a lot of memory!
#'
#' @param tree1 A phylogenetic tree (class \code{phylo}) or vector of trees (an
#' object of class \code{multiPhylo}). See details
#' @param tree2 A phylogenetic tree.
#' @param normalize compute normalized RF-distance, see details.
#' @param check.labels compares labels of the trees.
#' @param rooted take bipartitions for rooted trees into account, default is
#' unrooting the trees.
#' @param use.weight use edge.length argument or just count number of edges on
#' the path (default)
#' @return \code{treedist} returns a vector containing the following tree
#' distance methods \item{symmetric.difference}{symmetric.difference or
#' Robinson-Foulds distance}
#' \item{branch.score.difference}{branch.score.difference}
#' \item{path.difference}{path.difference}
#' \item{weighted.path.difference}{weighted.path.difference}
#' @author Klaus P. Schliep \email{klaus.schliep@@gmail.com},
#' Leonardo de Oliveira Martins
#' @seealso \code{\link[ape]{dist.topo}}, \code{\link{nni}},
#' \code{\link{superTree}}, \code{\link{mast}}
#' @references de Oliveira Martins L., Leal E., Kishino H. (2008)
#' \emph{Phylogenetic Detection of Recombination with a Bayesian Prior on the
#' Distance between Trees}. PLoS ONE \bold{3(7)}. e2651. doi:
#' 10.1371/journal.pone.0002651
#'
#' de Oliveira Martins L., Mallo D., Posada D. (2016) \emph{A Bayesian
#' Supertree Model for Genome-Wide Species Tree Reconstruction}. Syst. Biol.
#' \bold{65(3)}: 397-416, doi:10.1093/sysbio/syu082
#'
#' Steel M. A. and Penny P. (1993) \emph{Distributions of tree comparison
#' metrics - some new results}, Syst. Biol., \bold{42(2)}, 126--141
#'
#' Kuhner, M. K. and Felsenstein, J. (1994) \emph{A simulation comparison of
#' phylogeny algorithms under equal and unequal evolutionary rates}, Molecular
#' Biology and Evolution, \bold{11(3)}, 459--468
#'
#' D.F. Robinson and L.R. Foulds (1981) \emph{Comparison of phylogenetic
#' trees}, Mathematical Biosciences, \bold{53(1)}, 131--147
#'
#' D.F. Robinson and L.R. Foulds (1979) Comparison of weighted labelled trees.
#' In Horadam, A. F. and Wallis, W. D. (Eds.), \emph{Combinatorial Mathematics
#' VI: Proceedings of the Sixth Australian Conference on Combinatorial
#' Mathematics, Armidale, Australia}, 119--126
#' @keywords classif
#' @importFrom fastmatch fmatch
#' @examples
#'
#' tree1 <- rtree(100, rooted=FALSE)
#' tree2 <- rSPR(tree1, 3)
#' RF.dist(tree1, tree2)
#' treedist(tree1, tree2)
#' sprdist(tree1, tree2)
#' trees <- rSPR(tree1, 1:5)
#' SPR.dist(tree1, trees)
#'
#' @rdname treedist
#' @export treedist
treedist <- function(tree1, tree2, check.labels = TRUE) {
  if (has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if (has.singles(tree2)) tree2 <- collapse.singles(tree2)

  tree1 <- unroot(tree1)
  tree2 <- unroot(tree2)
  if (check.labels) tree2 <- checkLabels(tree2, tree1$tip.label)
#  if (check.labels) {
#    ind <- match(tree1$tip.label, tree2$tip.label)
#    if (any(is.na(ind)) | length(tree1$tip.label) !=
#      length(tree2$tip.label))
#      stop("trees have different labels")
#    tree2$tip.label <- tree2$tip.label[ind]
#    ind2 <- match(seq_along(ind), tree2$edge[, 2])
#    tree2$edge[ind2, 2] <- order(ind)
#  }

  tree1 <- reorder(tree1, "postorder")
  tree2 <- reorder(tree2, "postorder")

  symmetric.difference <- NULL
  branch.score.difference <- NULL
  path.difference <- NULL
  quadratic.path.difference <- NULL
  if (!is.binary(tree1) | !is.binary(tree2)) message("Trees are not binary!")

  bp1 <- bip(tree1)
  bp2 <- bip(tree2)
  bp1 <- SHORTwise(bp1)
  bp2 <- SHORTwise(bp2)
  bp1 <- sapply(bp1, paste, collapse = "_")
  bp2 <- sapply(bp2, paste, collapse = "_")

  l <- length(tree1$tip.label)

  if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {
    dv1 <- coph(tree1)
    dv2 <- coph(tree2)
    quadratic.path.difference <- sqrt(sum( (dv1 - dv2)^2))
  }

  RF <- sum(match(bp1, bp2, nomatch = 0L) == 0L) +
    sum(match(bp2, bp1, nomatch = 0L) == 0L)

  symmetric.difference <- RF # 2 * (p - sum(r1))
  if (!is.null(tree1$edge.length) & !is.null(tree2$edge.length)) {
    w1 <- numeric(max(tree1$edge))
    w2 <- numeric(max(tree2$edge))
    w1[tree1$edge[, 2]] <- tree1$edge.length
    w2[tree2$edge[, 2]] <- tree2$edge.length

#    v1 <- tree1$edge.length
#    v2 <- tree2$edge.length

    ind3 <- match(bp1, bp2, nomatch = 0L)
    ind4 <- ind3[ind3 > 0]
    ind3 <- which(ind3 > 0)

    s1 <- sum( (w1[ind3] - w2[ind4])^2)

    s2 <- sum(w1[-ind3]^2)
    s3 <- sum(w2[-ind4]^2)
    branch.score.difference <- sqrt(s1 + s2 + s3)
  }

  tree1$edge.length <- rep(1, nrow(tree1$edge))
  tree2$edge.length <- rep(1, nrow(tree2$edge))

  dt1 <- coph(tree1)
  dt2 <- coph(tree2)
  path.difference <- sqrt(sum( (dt1 - dt2)^2))

  result <- c(symmetric.difference = symmetric.difference,
    branch.score.difference = branch.score.difference,
    path.difference = path.difference,
    quadratic.path.difference = quadratic.path.difference)
  result
}




# leomrtns addition
#' @rdname treedist
#' @export
sprdist <- function(tree1, tree2) {
  if (has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if (has.singles(tree2)) tree2 <- collapse.singles(tree2)
  tree1 <- unroot(tree1)
  tree2 <- unroot(tree2)
  lt1 <- length(tree1$tip.label)
  lt2 <- length(tree2$tip.label)
  # checking labels is obligatory for spr (user should prune one of them
  # beforehand?)
  ind <- match(tree1$tip.label, tree2$tip.label)
  if (any(is.na(ind)) | lt1 != lt2) stop("trees have different labels")
  tree2$tip.label <- tree2$tip.label[ind]
  ind2 <- match(seq_along(ind), tree2$edge[, 2])
  tree2$edge[ind2, 2] <- order(ind)
  # same as in original treedist (will create list of strings with shorter
  # side of splits)
  tree1 <- reorder(tree1, "postorder")
  tree2 <- reorder(tree2, "postorder")
  if (!is.binary(tree1) | !is.binary(tree2)) message("Trees are not binary!")
  # possibly replace bip with bipart
  bp1 <- bip(tree1)
  bp1 <- SHORTwise(bp1)
  bp2 <- bip(tree2)
  bp2 <- SHORTwise(bp2)

  bp1 <- bp1[ lengths(bp1) > 1 ] # only internal nodes
  bp2 <- bp2[ lengths(bp2) > 1 ]
  if (length(bp1) != length(bp2))
    stop ("number of bipartitions given to C_sprdist are not the same")
  # OBS: SPR distance works w/ incompatible splits only, but it needs common
  # cherries (to replace by single leaf)
  spr <- .Call("C_sprdist", bp1, bp2, lt1)
  tmp <- .Call("C_sprdist", bp2, bp1, lt1)[1]
  spr[1] <- min(spr[1], tmp)
  names(spr) <- c("spr", "spr_extra", "rf", "hdist")
  spr
}


SPR1 <- function(trees) {
  trees <- .compressTipLabel(trees)
  trees <- .uncompressTipLabel(trees)
  trees <- lapply(trees, unroot)
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  trees <- lapply(trees, reorder, "postorder")

  nTips <- length(trees[[1]]$tip.label)

  fun <- function(x) {
    bp <- bipart(x)
    bp <- SHORTwise(bp)
    bp <- bp[ lengths(bp) > 1 ]
    bp
  }

  BP <- lapply(trees, fun)
  k <- 1
  l <- length(trees)
  SPR <- numeric( (l * (l - 1)) / 2)
  for (i in 1:(l - 1)) {
    bp <- BP[[i]]
    for (j in (i + 1):l) {
      SPR[k] <-  min( .Call("C_sprdist", bp, BP[[j]], nTips)[1],
                      .Call("C_sprdist", BP[[j]], bp, nTips)[1])
      k <- k + 1
    }
  }
  attr(SPR, "Size") <- l
  if (!is.null(names(trees))) attr(SPR, "Labels") <- names(trees)
  attr(SPR, "Diag") <- FALSE
  attr(SPR, "Upper") <- FALSE
  class(SPR) <- "dist"
  return(SPR)
}


SPR2 <- function(tree, trees) {
  trees <- .compressTipLabel(trees)
  tree <- checkLabels(tree, attr(trees, "TipLabel"))
  trees <- .uncompressTipLabel(trees)
  if (any(is.rooted(trees))) {
    trees <- unroot(trees)
  }
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  trees <- lapply(trees, reorder, "postorder")
  tree <- unroot(tree)
  if (has.singles(tree)) tree <- collapse.singles(tree)
  nTips <- length(tree$tip.label)

  fun <- function(x) {
    bp <- bipart(x)
    bp <- SHORTwise(bp)
    bp <- bp[ lengths(bp) > 1 ]
    bp
  }

  bp <-  fun(tree)
  l <- length(trees)
  SPR <- numeric(l)
  for (i in 1:l) {
    bpi <- fun(trees[[i]])
    SPR[i] <- min(.Call("C_sprdist", bp, bpi, nTips)[1],
                  .Call("C_sprdist", bpi, bp, nTips)[1])
  }
  if (!is.null(names(trees))) names(SPR) <- names(trees)
  return(SPR)
}


#' @rdname treedist
#' @export
SPR.dist <- function(tree1, tree2 = NULL) {
  if (inherits(tree1, "multiPhylo") && is.null(tree2)) return(SPR1(tree1))
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo"))
    return(sprdist(tree1, tree2)[1])
  if (inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
    return(SPR2(tree1, tree2))
  if (inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
    return(SPR2(tree2, tree1))
  return(NULL)
}


wRF0 <- function(tree1, tree2, normalize = FALSE, check.labels = TRUE,
                 rooted = FALSE) {
  r1 <- is.rooted(tree1)
  r2 <- is.rooted(tree2)
  if (r1 != r2) {
    message("one tree is unrooted, unrooted both")
  }
  if (!rooted) {
    if (r1)
      tree1 <- unroot(tree1)
    if (r2)
      tree2 <- unroot(tree2)
  }
  if (!r1 | !r2) {
    if (r1)
      tree1 <- unroot(tree1)
    if (r2)
      tree2 <- unroot(tree2)
  }
  if (!is.binary(tree1) | !is.binary(tree2))
    message("Trees are not binary!")
  if (check.labels) tree2 <- checkLabels(tree2, tree1$tip.label)
  if (has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if (has.singles(tree2)) tree2 <- collapse.singles(tree2)
  bp1 <- bip(tree1)
  bp2 <- bip(tree2)
  if (!rooted) {
    bp1 <- SHORTwise(bp1)
    bp2 <- SHORTwise(bp2)
  }
  bp1 <- sapply(bp1, paste, collapse = "_")
  bp2 <- sapply(bp2, paste, collapse = "_")

  w1 <- numeric(max(tree1$edge))
  w2 <- numeric(max(tree2$edge))
  w1[tree1$edge[, 2]] <- tree1$edge.length
  w2[tree2$edge[, 2]] <- tree2$edge.length

  ind3 <- match(bp1, bp2, nomatch = 0L)
  ind4 <- ind3[ind3 > 0]
  ind3 <- which(ind3 > 0)

  s1 <- sum(abs(w1[ind3] - w2[ind4]))
  s2 <- sum(w1[-ind3])
  s3 <- sum(w2[-ind4])

  wRF <- s1 + s2 + s3
  if (normalize) wRF <- wRF / (sum(tree1$edge.length) + sum(tree2$edge.length))
  return(wRF)
}


wRF2 <- function(tree, trees, normalize = FALSE, check.labels = TRUE,
                 rooted = FALSE) {
  if (check.labels) {
    trees <- .compressTipLabel(trees)
    tree <- checkLabels(tree, attr(trees, "TipLabel"))
  }
  trees <- .uncompressTipLabel(trees)

  if (rooted & any(!is.rooted(trees))) {
    warning("some trees were rooted, unrooted all")
    rooted <- FALSE
  }

  if (!rooted) {
    if (any(is.rooted(trees))) {
      trees <- unroot(trees)
    }
  }

  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  if (has.singles(tree)) tree <- collapse.singles(tree)

  unclass(trees)

  nTips <- length(tree$tip.label)

  fun1 <- function(x) {
    w <- numeric(max(x$edge))
    w[x$edge[, 2]] <- x$edge.length
    w
  }
  W <- lapply(trees, fun1)

  fun2 <- function(x, nTips) {
    bp <- bip(x)
    bp <- SHORTwise(bp)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  fun3 <- function(x, nTips) {
    bp <- bip(x)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  if (rooted) BP <- lapply(trees, fun3, nTips)
  else BP <- lapply(trees, fun2, nTips)

  if (!rooted & is.rooted(tree)) tree <- unroot(tree)

  bp <- bip(tree)

  if (!rooted) bp <- SHORTwise(bp)
  bp <- sapply(bp, paste, collapse = "_")

  w <- numeric(max(tree$edge))
  w[tree$edge[, 2]] <- tree$edge.length

  l <- length(trees)
  wRF <- numeric(l)

  for (i in 1:l) {
    ind3 <- fmatch(BP[[i]], bp, nomatch = 0L)
    ind4 <- ind3[ind3 > 0]
    ind3 <- which(ind3 > 0)

    s1 <- sum(abs(W[[i]][ind3] - w[ind4]))
    s2 <- sum(W[[i]][-ind3])
    s3 <- sum(w[-ind4])
    wRF[i] <- (s1 + s2 + s3)
  }
  if (normalize) {
    sc <- sapply(trees, function(x) sum(x$edge.length)) + sum(tree$edge.length)
    wRF <- wRF / sc
  }
  wRF
}


wRF1 <- function(trees, normalize = FALSE, check.labels = TRUE,
                 rooted = FALSE) {
  if (check.labels) trees <- .compressTipLabel(trees)
  trees <- .uncompressTipLabel(trees)

  if (rooted & any(!is.rooted(trees))) {
    warning("some trees were rooted, unrooted all")
    rooted <- FALSE
  }
  if (!rooted) {
    if (any(is.rooted(trees))) {
      trees <- unroot(trees)
    }
  }
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  unclass(trees)

  nTips <- length(trees[[1]]$tip.label)
  fun1 <- function(x) {
    w <- numeric(max(x$edge))
    w[x$edge[, 2]] <- x$edge.length
    w
  }
  W <- lapply(trees, fun1)
  fun2 <- function(x, nTips) {
    bp <- bip(x)
    bp <- SHORTwise(bp)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  fun3 <- function(x, nTips) {
    bp <- bip(x)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  if (normalize) sc <- sapply(trees, function(x) sum(x$edge.length))
  if (rooted) BP <- lapply(trees, fun3, nTips)
  else BP <- lapply(trees, fun2, nTips)
  k <- 1
  l <- length(trees)
  wRF <- numeric( (l * (l - 1)) / 2)
  for (i in 1:(l - 1)) {
    bp <- BP[[i]]
    w <- W[[i]]
    for (j in (i + 1):l) {
      ind3 <- fmatch(BP[[j]], bp, nomatch = 0L)
      ind4 <- ind3[ind3 > 0]
      ind3 <- which(ind3 > 0)
      s1 <- sum(abs(W[[j]][ind3] - w[ind4]))
      s2 <- sum(W[[j]][-ind3])
      s3 <- sum(w[-ind4])
      wRF[k] <- (s1 + s2 + s3)
      if (normalize) wRF[k] <- wRF[k] / (sc[i] + sc[j])
      k <- k + 1
    }
  }
  attr(wRF, "Size") <- l
  if (!is.null(names(trees))) attr(wRF, "Labels") <- names(trees)
  attr(wRF, "Diag") <- FALSE
  attr(wRF, "Upper") <- FALSE
  class(wRF) <- "dist"
  return(wRF)
}


mRF2 <- function(tree, trees, normalize = FALSE, check.labels = TRUE,
                 rooted = FALSE) {
  if (!inherits(trees, "multiPhylo"))
    stop("Argument trees should be an object of class \"multiPhylo\"")
  if (!inherits(tree, "phylo"))
    stop("Argument tree should be an object of class \"phylo\"")
  trees <- .compressTipLabel(trees)
  tipLabel <- attr(trees, "TipLabel")
  if (check.labels) tree <- checkLabels(tree, tipLabel)
#  if (check.labels) {
#    ind <- match(tipLabel, tree$tip.label)
#    if (any(is.na(ind)) | length(tipLabel) != length(tree$tip.label))
#      stop("trees have different labels")
#    tree$tip.label <- tree$tip.label[ind]
#    ind2 <- match(seq_along(ind), tree$edge[, 2])
#    tree$edge[ind2, 2] <- order(ind)
#  }
  nTips <- length(tipLabel)
  l <- length(trees)
  RF <- numeric(l)
  trees <- .uncompressTipLabel(trees)
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  if (has.singles(tree)) tree <- collapse.singles(tree)

  if (!rooted & any(is.rooted(trees))) {
    warning("some trees were rooted, unrooted all")
    trees <- unroot(trees)
  }
  if (!rooted & is.rooted(tree)) tree <- unroot(tree)
  if (any(!is.binary(trees))) {
    message("Some trees are not binary. Result may not what you expect!")
  }
  tree <- reorder(tree, "postorder")
  trees <- reorder(trees, "postorder")
  xx <- lapply(trees, bipart)
  if (!rooted) xx <- lapply(xx, SHORTwise)
  xx <- lapply(xx, function(x) sapply(x, paste, collapse = "_"))
  yy <- bipart(tree)
  if (!rooted) yy <- SHORTwise(yy)
  yy <- sapply(yy, paste, collapse = "_")

  NnodeT <- Nnode(tree)
  Nnodes <- Nnode(trees)

  for (i in 1:l) {
    RF[i] <- Nnodes[i] + NnodeT - 2 *
      sum(fmatch(xx[[i]], yy, nomatch = 0L) > 0L)
    # RF[i] <- sum(match(xx[[i]], yy, nomatch=0L)==0L) +
    #   sum(match(yy, xx[[i]], nomatch=0L)==0L)
  }
  if (!is.null(names(trees))) names(RF) <- names(trees)
  if (!normalize) return(RF)
  else {
    sc <- Nnode(trees) + Nnode(tree) - 2
    return(RF / sc)
  }
}


mRF <- function(trees, normalize = FALSE, rooted = FALSE) {
  if (!inherits(trees, "multiPhylo"))
    stop("Argument trees should be an object of class \"multiPhylo\"")
  trees <- .compressTipLabel(trees)
  tipLabel <- attr(trees, "TipLabel")
  nTips <- length(tipLabel)
  l <- length(trees)
  RF <- numeric( (l * (l - 1)) / 2)

  if (rooted & any(!is.rooted(trees))) {
    warning("Some trees were rooted, unrooted all")
    rooted <- FALSE
  }
  if (!rooted) {
    if (any(is.rooted(trees))) {
      trees <- unroot(trees)
    }
  }

  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)

  #    n <- length(attr(trees, "TipLabel"))

  #    if (any(sapply(trees, is.rooted))) {
  #        warning("some trees were rooted, unrooted all")
  #        trees <- lapply(trees, unroot)
  #    }
  if (any(!is.binary(trees))) {
    message("Some trees are not binary. Result may not what you expect!")
  }
  #    trees <- reorder(trees, "postorder")
  #    trees <- lapply(trees, reorder, "postorder")
  Nnodes <- Nnode(trees)
  trees <- .uncompressTipLabel(trees)
  trees <- unclass(trees)

  xx <- lapply(trees, bipart)
  if (!rooted) xx <- lapply(xx, SHORTwise)
  xx <- lapply(xx, function(x) sapply(x, paste, collapse = "_"))
  # returns list of character vectors

  k <- 1
  for (i in 1:(l - 1)) {
    tmp <- xx[[i]]
    for (j in (i + 1):l) {
      RF[k] <- Nnodes[i] + Nnodes[j] - 2 *
        sum(fmatch(xx[[j]], tmp, nomatch = 0L) > 0L)
      #  RF[k] <- sum(match(xx[[j]], tmp, nomatch=0L)==0L) +
      #  sum(match(tmp, xx[[j]], nomatch=0L)==0L)
      if (normalize) RF[k] <- RF[k] / (Nnodes[i] + Nnodes[j] - 2)
      k <- k + 1
    }
  }
  attr(RF, "Size") <- l
  if (!is.null(names(trees))) attr(RF, "Labels") <- names(trees)
  attr(RF, "Diag") <- FALSE
  attr(RF, "Upper") <- FALSE
  class(RF) <- "dist"
  return(RF)
}


RF0 <- function(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE,
                rooted = FALSE) {
  if (has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if (has.singles(tree2)) tree2 <- collapse.singles(tree2)
  r1 <- is.rooted(tree1)
  r2 <- is.rooted(tree2)
  if (!rooted) {
    if (r1) {
      tree1 <- unroot(tree1)
      r1 <- FALSE
    }
    if (r2) {
      tree2 <- unroot(tree2)
      r2 <- FALSE
    }
  }
  else {
    if (r1 != r2) {
      message("one tree is unrooted, unrooted both")
      tree1 <- unroot(tree1)
      tree2 <- unroot(tree2)
      r1 <- r2 <- FALSE
    }
  }
  if (check.labels) tree2 <- checkLabels(tree2, tree1$tip.label)
  if (!is.binary(tree1) | !is.binary(tree2)) message("Trees are not binary!")
  bp1 <- bipart(tree1)
  bp2 <- bipart(tree2)
  nTips <- length(tree1$tip.label)
  if (!rooted) {
    bp1 <- SHORTwise(bp1)
    bp2 <- SHORTwise(bp2)
  }
  RF <- sum(match(bp1, bp2, nomatch = 0L) == 0L) +
    sum(match(bp2, bp1, nomatch = 0L) == 0L)
  if (normalize) RF <- RF / (Nnode(tree1) + Nnode(tree2) - 2)
  RF
}


#' @rdname treedist
#' @export
RF.dist <- function(tree1, tree2 = NULL, normalize = FALSE, check.labels = TRUE,
                    rooted = FALSE) {
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo"))
    return(RF0(tree1, tree2, normalize, check.labels, rooted))
  if (inherits(tree1, "multiPhylo") && is.null(tree2))
    return(mRF(tree1, normalize, rooted))
  if (inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
    return(mRF2(tree1, tree2, normalize, check.labels, rooted))
  if (inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
    return(mRF2(tree2, tree1, normalize, check.labels, rooted))
  else return(NULL)
}


#' @rdname treedist
#' @export
wRF.dist <- function(tree1, tree2 = NULL, normalize = FALSE,
                     check.labels = TRUE, rooted = FALSE) {
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo"))
    return(wRF0(tree1, tree2, normalize, check.labels, rooted))
  if (inherits(tree1, "multiPhylo") && is.null(tree2))
    return(wRF1(tree1, normalize, check.labels, rooted))
  if (inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
    return(wRF2(tree1, tree2, normalize, check.labels, rooted))
  if (inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
    return(wRF2(tree2, tree1, normalize, check.labels, rooted))
  else return(NULL)
}


kf0 <- function(tree1, tree2, check.labels = TRUE, rooted = FALSE) {
  if (check.labels) tree2 <- checkLabels(tree2, tree1$tip.label)
  if (has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if (has.singles(tree2)) tree2 <- collapse.singles(tree2)
  r1 <- is.rooted(tree1)
  r2 <- is.rooted(tree2)
  if (!rooted) {
    if (r1) tree1 <- unroot(tree1)
    if (r2) tree2 <- unroot(tree2)
  }
  else {
    if (r1 != r2) {
      message("one tree is unrooted, unrooted both")
      tree1 <- unroot(tree1)
      tree2 <- unroot(tree2)
      r1 <- r2 <- FALSE
    }
  }

  bp1 <- bip(tree1)
  bp2 <- bip(tree2)

  if (!rooted) {
    bp1 <- SHORTwise(bp1)
    bp2 <- SHORTwise(bp2)
  }
  bp1 <- sapply(bp1, paste, collapse = "_")
  bp2 <- sapply(bp2, paste, collapse = "_")

  w1 <- numeric(max(tree1$edge))
  w2 <- numeric(max(tree2$edge))
  w1[tree1$edge[, 2]] <- tree1$edge.length
  w2[tree2$edge[, 2]] <- tree2$edge.length

  ind3 <- match(bp1, bp2, nomatch = 0L)
  ind4 <- ind3[ind3 > 0]
  ind3 <- which(ind3 > 0)

  s1 <- sum( (w1[ind3] - w2[ind4])^2)
  s2 <- sum(w1[-ind3]^2)
  s3 <- sum(w2[-ind4]^2)
  branch.score.difference <- sqrt(s1 + s2 + s3)
  branch.score.difference
}


kf1 <- function(tree, trees, check.labels = TRUE, rooted = FALSE) {
  if (check.labels) {
    trees <- .compressTipLabel(trees)
    tree <- checkLabels(tree, attr(trees, "TipLabel"))
  }
  trees <- .uncompressTipLabel(trees)
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)
  if (has.singles(tree)) tree <- collapse.singles(tree)

  if (rooted & any(!is.rooted(trees))) {
    warning("some trees were rooted, unrooted all")
    rooted <- FALSE
  }
  if (!rooted) {
    if (any(is.rooted(trees))) {
      trees <- unroot(trees)
    }
  }

  unclass(trees)

  nTips <- length(tree$tip.label)

  fun1 <- function(x) {
    w <- numeric(max(x$edge))
    w[x$edge[, 2]] <- x$edge.length
    w
  }
  W <- lapply(trees, fun1)

  fun2 <- function(x, nTips) {
    bp <- bip(x)
    bp <- SHORTwise(bp)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  fun3 <- function(x, nTips) {
    bp <- bip(x)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  if (rooted) BP <- lapply(trees, fun3, nTips)
  else BP <- lapply(trees, fun2, nTips)

  if (!rooted & is.rooted(tree)) tree <- unroot(tree)
  bp <- bip(tree)
  if (!rooted) bp <- SHORTwise(bp)
  bp <- sapply(bp, paste, collapse = "_")

  w <- numeric(max(tree$edge))
  w[tree$edge[, 2]] <- tree$edge.length

  l <- length(trees)
  branch.score.difference <- numeric(l)

  for (i in 1:l) {
    ind3 <- fmatch(BP[[i]], bp, nomatch = 0L)
    ind4 <- ind3[ind3 > 0]
    ind3 <- which(ind3 > 0)

    s1 <- sum( (W[[i]][ind3] - w[ind4])^2)
    s2 <- sum(W[[i]][-ind3]^2)
    s3 <- sum(w[-ind4]^2)
    branch.score.difference[i] <- sqrt(s1 + s2 + s3)
  }
  branch.score.difference
}


kf2 <- function(trees, check.labels = TRUE, rooted = FALSE) {
  if (check.labels) trees <- .compressTipLabel(trees)
  trees <- .uncompressTipLabel(trees)
  if (any(has.singles(trees))) trees <- lapply(trees, collapse.singles)

  nTips <- length(trees[[1]]$tip.label)
  if (rooted & any(!is.rooted(trees))) {
    warning("some trees were rooted, unrooted all")
    rooted <- FALSE
  }
  if (!rooted & any(is.rooted(trees))) {
    trees <- unroot(trees)
  }

  unclass(trees)
  fun1 <- function(x) {
    w <- numeric(max(x$edge))
    w[x$edge[, 2]] <- x$edge.length
    w
  }
  W <- lapply(trees, fun1)


  fun2 <- function(x, nTips) {
    bp <- bip(x)
    bp <- SHORTwise(bp)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  fun3 <- function(x, nTips) {
    bp <- bip(x)
    bp <- sapply(bp, paste, collapse = "_")
    bp
  }
  if (rooted) BP <- lapply(trees, fun3, nTips)
  else BP <- lapply(trees, fun2, nTips)

  k <- 1
  l <- length(trees)
  KF <- numeric( (l * (l - 1)) / 2)
  for (i in 1:(l - 1)) {
    bp <- BP[[i]]
    w <- W[[i]]
    for (j in (i + 1):l) {
      ind3 <- fmatch(BP[[j]], bp, nomatch = 0L)
      ind4 <- ind3[ind3 > 0]
      ind3 <- which(ind3 > 0)
      s1 <- sum( (W[[j]][ind3] - w[ind4])^2)
      s2 <- sum(W[[j]][-ind3]^2)
      s3 <- sum(w[-ind4]^2)
      KF[k] <- sqrt(s1 + s2 + s3)
      k <- k + 1
    }
  }
  attr(KF, "Size") <- l
  if (!is.null(names(trees))) attr(KF, "Labels") <- names(trees)
  attr(KF, "Diag") <- FALSE
  attr(KF, "Upper") <- FALSE
  class(KF) <- "dist"
  return(KF)
}


#' @rdname treedist
#' @export
KF.dist <- function(tree1, tree2 = NULL, check.labels = TRUE, rooted = FALSE) {
  if (inherits(tree1, "multiPhylo") && is.null(tree2))
    return(kf2(tree1, rooted = rooted))
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo"))
    return(kf0(tree1, tree2, check.labels, rooted))
  if (inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
    return(kf1(tree1, tree2, check.labels, rooted))
  if (inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
    return(kf1(tree2, tree1, check.labels, rooted))
  return(NULL)
}


#' @rdname treedist
#' @export
path.dist <- function(tree1, tree2 = NULL, check.labels = TRUE,
                      use.weight = FALSE) {
  if (inherits(tree1, "phylo") && inherits(tree2, "phylo"))
    return(pd0(tree1, tree2, check.labels, !use.weight))
  if (inherits(tree1, "phylo") && inherits(tree2, "multiPhylo"))
    return(pd1(tree1, tree2, check.labels, !use.weight))
  if (inherits(tree2, "phylo") && inherits(tree1, "multiPhylo"))
    return(pd1(tree2, tree1, check.labels, !use.weight))
  if (inherits(tree1, "multiPhylo") && is.null(tree2))
    return(pd2(tree1, check.labels, !use.weight))
  else return(NULL)
}


pd0 <- function(tree1, tree2, check.labels = TRUE, path = TRUE) {
  if (check.labels) tree2 <- checkLabels(tree2, tree1$tip.label)
  if (path) {
    tree1 <- unroot(tree1)
    tree2 <- unroot(tree2)
  }
  dt1 <- coph(tree1, path)
  dt2 <- coph(tree2, path)
  sqrt(sum( (dt1 - dt2)^2))
}


pd1 <- function(tree, trees, check.labels = TRUE, path = TRUE) {
  if (check.labels) {
    trees <- .compressTipLabel(trees)
    tree <- checkLabels(tree, attr(trees, "TipLabel"))
  }
  trees <- .uncompressTipLabel(trees)
  if (path) {
    trees <- unroot(trees)
    tree <- unroot(tree)
  }
  trees <- reorder(trees, "postorder")
  unclass(trees)
  l <- length(trees)
  dt <- coph(tree, path)
  res <- numeric(l)
  for (i in 1:l) {
    dt2 <- coph(trees[[i]], path)
    res[i] <- sqrt(sum( (dt - dt2)^2))
  }
  res
}

pd2 <- function(trees, check.labels = TRUE, path = TRUE) {
  if (check.labels) trees <- .compressTipLabel(trees)
  trees <- .uncompressTipLabel(trees)
  if (path) trees <- unroot(trees)
  trees <- reorder(trees, "postorder")
  l <- length(trees)
  unclass(trees)
  CM <- lapply(trees, coph, path)
  k <- 1
  PD <- numeric( (l * (l - 1)) / 2)
  for (i in 1:(l - 1)) {
    for (j in (i + 1):l) {
      PD[k] <- sqrt(sum( (CM[[i]] - CM[[j]])^2))
      k <- k + 1
    }
  }
  attr(PD, "Size") <- l
  if (!is.null(names(trees))) attr(PD, "Labels") <- names(trees)
  attr(PD, "Diag") <- FALSE
  attr(PD, "Upper") <- FALSE
  class(PD) <- "dist"
  return(PD)
}
