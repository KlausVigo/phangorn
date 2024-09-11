#' Parsimony tree.
#'
#' \code{pratchet} implements the parsimony ratchet (Nixon, 1999) and is the
#' preferred way to search for the best parsimony tree. For small number of taxa
#' the function \code{\link{bab}} can be used to compute all most parsimonious
#' trees.
#'
#' \code{parsimony} returns the parsimony score of a tree using either the
#' sankoff or the fitch algorithm.
#' \code{optim.parsimony} optimizes the topology using either Nearest Neighbor
#' Interchange (NNI) rearrangements or sub tree pruning and regrafting (SPR) and
#' is used inside \code{pratchet}. \code{random.addition} can be used to produce
#' starting trees and is an option for the argument perturbation in
#' \code{pratchet}.
#'
#' The "SPR" rearrangements are so far only available for the "fitch" method,
#' "sankoff" only uses "NNI". The "fitch" algorithm only works correct for
#' binary trees.
#'
#' @aliases parsimony
#' @param data A object of class phyDat containing sequences.
#' @param tree tree to start the nni search from.
#' @param method one of 'fitch' or 'sankoff'.
#' @param cost A cost matrix for the transitions between two states.
#' @param site return either 'pscore' or 'site' wise parsimony scores.
#' @param trace defines how much information is printed during optimization.
#' @param rearrangements SPR or NNI rearrangements.
#' @param start a starting tree can be supplied.
#' @param maxit maximum number of iterations in the ratchet.
#' @param minit minimum number of iterations in the ratchet.
#' @param k number of rounds ratchet is stopped, when there is no improvement.
#' @param all return all equally good trees or just one of them.
#' @param perturbation whether to use "ratchet", "random_addition" or
#' "stochastic" (nni) for shuffling the tree.

#' @param ... Further arguments passed to or from other methods (e.g.
#' model="sankoff" and cost matrix).
#' @return \code{parsimony} returns the maximum parsimony score (pscore).
#' \code{optim.parsimony} returns a tree after NNI rearrangements.
#' \code{pratchet} returns a tree or list of trees containing the best tree(s)
#' found during the search.  \code{acctran} returns a tree with edge length
#' according to the ACCTRAN criterion.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{bab}}, \code{\link{CI}}, \code{\link{RI}},
#' \code{\link{ancestral.pml}}, \code{\link{nni}}, \code{\link{NJ}},
#' \code{\link{pml}}, \code{\link{getClans}}, \code{\link{ancestral.pars}},
#' \code{\link{bootstrap.pml}}
#' @references Felsenstein, J. (2004). \emph{Inferring Phylogenies}. Sinauer
#' Associates, Sunderland.
#'
#' Nixon, K. (1999) The Parsimony Ratchet, a New Method for Rapid Parsimony
#' Analysis. \emph{Cladistics} \bold{15}, 407-414
#' @keywords cluster
#' @examples
#'
#' set.seed(3)
#' data(Laurasiatherian)
#' dm <- dist.hamming(Laurasiatherian)
#' tree <- NJ(dm)
#' parsimony(tree, Laurasiatherian)
#' treeRA <- random.addition(Laurasiatherian)
#' treeSPR <- optim.parsimony(tree, Laurasiatherian)
#'
#' # lower number of iterations for the example (to run less than 5 seconds),
#' # keep default values (maxit, minit, k) or increase them for real life
#' # analyses.
#' treeRatchet <- pratchet(Laurasiatherian, start=tree, maxit=100,
#'                         minit=5, k=5, trace=0)
#' # assign edge length (number of substitutions)
#' treeRatchet <- acctran(treeRatchet, Laurasiatherian)
#' # remove edges of length 0
#' treeRatchet <- di2multi(treeRatchet)
#'
#' plot(midpoint(treeRatchet))
#' add.scale.bar(0,0, length=100)
#'
#' parsimony(c(tree,treeSPR, treeRatchet), Laurasiatherian)
#'
#' @rdname parsimony
#' @export
## parsimony <- function(tree, data, cost=NULL, method = NULL)
parsimony <- function(tree, data, method = "fitch", cost=NULL, site = "pscore"){
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  method <- match.arg(method, c("fitch", "sankoff"))
  if(!any(is.binary(tree)) || !is.null(cost)) method <- "sankoff"
  if (method == "sankoff") result <- sankoff(tree, data, cost=cost, site = site)
  if (method == "fitch") result <- fitch(tree, data, site = site)
  result
}


compressSites <- function(data) {
  attrData <- attributes(data)
  lev <- attr(data, "levels")
  LEV <- attr(data, "allLevels")
  l <- length(lev)
  nr <- attr(data, "nr")
  nc <- length(data)
  data <- unlist(data, FALSE, FALSE)
  attr(data, "dim") <- c(nr, nc)
  uni <- match(lev, LEV)
  fun <- function(x, uni) {
    u <- unique.default(x)
    res <- if (any(is.na(match(u, uni)))) return(x)
    match(x, u)
  }
  data <- apply(data, 1, fun, uni)
  index <- grp_duplicated(data, MARGIN=2L)
  pos <- which(!duplicated(index))
  ntaxa <- nrow(data)
  res <- vector("list", ntaxa)
  for(i in seq_len(ntaxa)) res[[i]] <- data[i, pos]
  attrData$weight <- as.vector(tapply(attrData$weight, index, sum))
  attrData$index <- NULL
  attrData$nr <- length(attrData$weight)
  attrData$compressed <- TRUE
  attributes(res) <- attrData
  res
}


parsinfo <- function(x, exact=TRUE) {
  nstates <- attr(x, "nc")
  up <- upperBound(x)
  eps <- 1e-8
  low <- up
  low[up > (nstates - eps)] <- nstates - 1
  if(exact){
    ind <- which( (up > (1+eps))  & (up < (nstates-eps)) )
    if(length(ind)>0) low[ind] <- lowerBound(getRows( x, ind ))
    ind <- which(low == up)
  }
  else ind <- which(up < (1+eps) )
  cbind(ind, low[ind])
}


# greedy algorithm of Maximum Set Packing (MSP) problem
# (should work in most instances)
lowerBound <- function(x, cost = NULL) {
  nc <- attr(x, "nc")
  nr <- attr(x, "nr")
  contrast <- attr(x, "contrast")
  rownames(contrast) <- attr(x, "allLevels")
  colnames(contrast) <- attr(x, "levels")
  nmax <- nrow(contrast)
  z <- matrix(unlist(x, FALSE, FALSE), length(x), length(attr(x, "weight")),
              byrow = TRUE)
  states <- apply(z, 2, unique.default, nmax = nmax)
  if(inherits(states, "matrix"))states <- asplit(states, 2)
  singles <- which(rowSums(contrast) == 1)
  noinfo <- which(rowSums(contrast) == nc)
  ambiguous <- which( (rowSums(contrast) > 1) & (rowSums(contrast) < nc))

  fun <- function(states, contrast, singles, noinfo, ambiguous) {
    if (length(states) == 1) return(0)
    states <- setdiff(states, noinfo) # get rid of "-", "?" in DNA
    if ( (length(states) == 0) | (length(states) == 1)) return(0)
    if (any(states %in% ambiguous)) {
      n <- 0L
      contrast <- contrast[states, , drop = FALSE]
      while (nrow(contrast) > 0) {
        m <- which.max(colSums(contrast))
        contrast <- contrast[contrast[, m] == 0, , drop = FALSE]
        n <- n + 1L
      }
      return(n - 1L)
    }
    else return(length(states) - 1L)
  }
  res <- sapply(states, fun, contrast, singles, noinfo, ambiguous)
  res
}


upperBound <- function(x, cost = NULL) {
  tree <- stree(length(x), tip.label = names(x))
  if (is.null(cost)) cost <- 1 - diag(attr(x, "nc"))
  sankoff(tree, x, cost = cost, site = "site")
}



#' Consistency Index and Retention Index
#'
#' \code{CI} and \code{RI} compute the Consistency Index (CI) and Retention
#' Index (RI).
#'
#' @details The Consistency Index is defined as minimum number of changes
#' divided by the number of changes required on the tree (parsimony score). The
#' Consistency Index is equal to one if there is no homoplasy.
#' And the Retention Index is defined as
#' \deqn{RI = \frac{MaxChanges - ObsChanges}{MaxChanges - MinChanges}}{RI = (MaxChanges - ObsChanges) / (MaxChanges - MinChanges)}
#'
#' @param data A object of class phyDat containing sequences.
#' @param tree a phylogenetic tree, i.e. an object of class \code{phylo}.
#' @param cost A cost matrix for the transitions between two states.
#' @param sitewise return CI/RI for alignment or sitewise
#' @returns a scalar or vector with the CI/RI vector.
#' @seealso \code{\link{parsimony}}, \code{\link{pratchet}},
#' \code{\link{fitch}}, \code{\link{sankoff}}, \code{\link{bab}},
#' \code{\link{ancestral.pars}}
#' @examples
#' example(as.phylo.formula)
#' lab <- tr$tip.label
#' lab[79] <- "Herpestes fuscus"
#' tr$tip.label <- abbreviateGenus(lab)
#' X <- matrix(0, 112, 3, dimnames = list(tr$tip.label, c("Canis", "Panthera",
#'             "Canis_Panthera")))
#' desc_canis <- Descendants(tr, "Canis")[[1]]
#' desc_panthera <- Descendants(tr, "Panthera")[[1]]
#' X[desc_canis, c(1,3)] <- 1
#' X[desc_panthera, c(2,3)] <- 1
#' col <- rep("black", 112)
#' col[desc_panthera] <- "red"
#' col[desc_canis] <- "blue"
#' X <- phyDat(X, "USER", levels=c(0,1))
#' plot(tr, "f", tip.color=col)
#' # The first two sites are homoplase free!
#' CI(tr, X, sitewise=TRUE)
#' RI(tr, X, sitewise=TRUE)
#'
#' @rdname CI
#' @export
CI <- function(tree, data, cost = NULL, sitewise = FALSE) {
  data <- subset(data, tree$tip.label)
  pscore <- sankoff(tree, data, cost, ifelse(sitewise, "site", "pscore"))
  weight <- attr(data, "weight")
  m <- lowerBound(data, cost = cost)
  if (sitewise) {
    return( (m / pscore)[attr(data, "index")])
  }
  sum(m * weight) / pscore
}


#' @rdname CI
#' @export
RI <- function(tree, data, cost = NULL, sitewise = FALSE) {
  data <- subset(data, tree$tip.label)
  pscore <- sankoff(tree, data, cost, ifelse(sitewise, "site", "pscore"))
  weight <- attr(data, "weight")
  m <- lowerBound(data, cost = cost)
  g <- upperBound(data, cost = cost)
  if (sitewise) {
    res <- (g - pscore) / (g - m)
    return(res[attr(data, "index")])
  }
  m <- sum(m * weight)
  g <- sum(g * weight)
  (g - pscore) / (g - m)
}


old2new.phyDat <- function(obj) {
  att <- attributes(obj)
  l <- length(obj)
  contrast <- attr(obj, "contrast")
  nr <- attr(obj, "nr")
  X <- matrix(rep(rowSums(contrast), each = nr), nrow = nr)
  for (i in 1:l) obj[[i]][obj[[i]] > 0] <- 1
  res <- vector("list", l)
  contrast[contrast == 0] <- 1e6
  for (i in 1:l) {
    tmp <-  tcrossprod(obj[[i]], contrast) - X
    res[[i]] <- unlist(apply(tmp, 1, function(x) which(x < 1e-6)[1]))
  }
  attributes(res) <- att
  res
}



new2old.phyDat <- function(data) {
  contrast <- attr(data, "contrast")
  for (i in seq_along(data)) data[[i]] <- contrast[data[[i]], , drop = FALSE]
  data
}


indexNNI <- function(tree) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]

  ind <- which(child %in% parent)
  Nnode <- tree$Nnode
  edgeMatrix <- matrix(0, (Nnode - 1), 5)

  pvector <- integer(max(parent))
  pvector[child] <- parent
  tips  <- !logical(max(parent))
  tips[parent] <-  FALSE
  # wahrscheinlich schneller: cvector <- allChildren(tree)
  cvector <- vector("list", max(parent))
  for (i in seq_along(parent)) cvector[[parent[i]]] <- c(cvector[[parent[i]]],
      child[i])
  k <- 0
  for (i in ind) {
    p1 <- parent[i]
    p2 <- child[i]
    e34 <- cvector[[p2]]
    ind1 <- cvector[[p1]]
    e12 <- ind1[ind1 != p2]
    if (pvector[p1]) e12 <- c(p1, e12)
    edgeMatrix[k + 1, ] <- c(e12, e34, p2)
    k <- k + 1
  }
  # vielleicht raus
  attr(edgeMatrix, "root") <- cvector[[min(parent)]]
  edgeMatrix
}


#' @rdname parsimony
#' @export
optim.parsimony <- function(tree, data, method = "fitch", cost = NULL,
                            trace = 1, rearrangements = "SPR", ...) {
  if (method == "fitch") result <- optim.fitch(tree = tree, data = data,
                      trace = trace, rearrangements = rearrangements, ...)
  if (method == "sankoff") result <- optim.sankoff(tree = tree, data = data,
      cost = cost, trace = trace, ...)
  result
}


#' @rdname parsimony
#' @export
pratchet <- function(data, start = NULL, method = "fitch", maxit = 1000,
                     minit = 100, k = 10, trace = 1, all = FALSE,
                     rearrangements = "SPR", perturbation = "ratchet", ...) {
  if(inherits(data, "DNAbin") || inherits(data, "AAbin"))
    data <- as.phyDat(data)
  printevery <- 10L
  eps <- 1e-08
  trace <- trace - 1
  ref <- names(data)
  start_trees <- vector("list", maxit)
  search_trees <- vector("list", maxit)
  tree <- NULL
  mp <- Inf
  # TODO use rooted trees if cost is not symmetric
  ROOTED <- FALSE
  weight <- attr(data, "weight")
  v <- rep(seq_along(weight), weight)
  w <- logical(length(weight))
  # remove parsimony uniformative sites or duplicates
  # check for symmetric or
  attr(data, "informative") <- NULL
  if(method=="fitch") data <- removeParsimonyUninfomativeSites(data,
                                                               recursive=TRUE)
  else data <- unique(data)
  if(!is.null(attr(data, "informative"))) w[attr(data, "informative")] <- TRUE
  else w[] <- TRUE

  star_tree <- ifelse(attr(data, "nr") == 0, TRUE, FALSE)
  add_taxa <- ifelse(is.null(attr(data, "duplicated")), FALSE, TRUE)
  nTips <- length(data)
  # check for trivial trees
  if (nTips < (3L + !ROOTED)  || star_tree) {
    nam <- names(data)
    if (star_tree) tree <- stree(length(nam), tip.label = nam)
    else tree <- stree(nTips, tip.label = nam)
    if(add_taxa) tree <- addTaxa(tree, attr(data, "duplicated"))
    if(!ROOTED) tree <- unroot(tree)
    tree <- relabel(tree, ref)
    return(tree)
  }

    if(is.null(start)) start <- optim.parsimony(random.addition(data),
                                   data, trace = trace-1, method = method,
                                   rearrangements = rearrangements, ...)
    tree <- start
    label <- intersect(tree$tip.label, names(data))
    if (!is.binary(tree)){
      tree <- multi2di(tree)
      if(method=="fitch") tree <- unroot(tree)
    }
    data <- subset(data, label)
    tree <- keep.tip(tree, label)
    attr(tree, "pscore") <- parsimony(tree, data, method = method, ...)
    mp <- attr(tree, "pscore")
    if (trace >= 0)
      cat("Parsimony score of initial tree:", attr(tree, "pscore"), "\n")
  FUN <- function(data, tree, method, rearrangements, ...)
    optim.parsimony(tree, data = data, method = method,
                    rearrangements = rearrangements, ...)
  result <- tree
  if(!is.null(attr(data, "duplicated"))){
    result <- addTaxa(result, attr(data, "duplicated"))
  }
  result <- relabel(result, ref)
#  if (trace > 1) cat("optimize topology (NNI): ", pscore, "-->", psc, "\n")
  hr <- hash(result)
  on.exit({
    if (!all && inherits(result, "multiPhylo")) result <- result[[1]]
    if (length(result) == 1) result <- result[[1]]
    env <- new.env()
    start_trees <- start_trees[seq_len(i)]
    search_trees <- search_trees[seq_len(i)]
    class(start_trees) <- "multiPhylo"
    class(search_trees) <- "multiPhylo"
    start_trees <- .compressTipLabel(start_trees)
    search_trees <- .compressTipLabel(search_trees)
    assign("start_trees", start_trees, envir=env)
    assign("search_trees", search_trees, envir=env)
    if(perturbation == "ratchet" &&  all(Ntip(trees) > 3)) {
      spl <- as.splits(start_trees)
      result <- addConfidences(result, spl)
      if (inherits(result, "multiPhylo")) result <- .compressTipLabel(result)
    }
    # for ratchet assign bs values
    attr(result, "env") <- env
    return(result)
  })
  kmax <- 1
  nTips <- length(tree$tip.label)
  for (i in seq_len(maxit)) {
    if (perturbation == "ratchet") {
      # sample and subset more efficient than in bootstrap.phyDat
      bsw <- tabulate(sample(v, replace = TRUE), length(weight))[w]
      bs_ind <- which(bsw > 0)
      bs_data <- getRows(data, bs_ind)
      attr(bs_data, "weight") <- bsw[bs_ind]
      if(length(bs_ind) > 0){
        # p_trees <- random.addition(bs_data)  # 3 * ??
        p_trees <- optim.parsimony(tree, bs_data,
          trace = trace, method = method, rearrangements = rearrangements, ...)
      }
      else p_trees <- stree(length(data), tip.label = names(data))
    }
    if (perturbation == "stochastic") p_trees <- rNNI(tree, floor(nTips / 2))
    if (perturbation == "random_addition") p_trees <- random.addition(data)
    trees <- optim.parsimony(p_trees, data, trace = trace, method = method,
                             rearrangements = rearrangements, ...)
    curr_tree <- trees
    if(!is.null(attr(data, "duplicated"))){
      p_trees <- addTaxa(p_trees, attr(data, "duplicated"))
      trees <- addTaxa(trees, attr(data, "duplicated"))
    }
    trees <- relabel(trees, ref)
    p_trees <- relabel(p_trees, ref)
    start_trees[[i]] <- p_trees
    search_trees[[i]] <- trees
    pscores <- attr(trees, "pscore")
    mp1 <- min(pscores)
    if ( (mp1 + eps) < mp) {
      kmax <- 1
      result <- trees
      tree <- curr_tree
      hr <- hash(trees)
      mp <- mp1
    }
    else{
      kmax <- kmax + 1
      if( all && (mp1 < (mp + eps))){ # && all(RF.dist(trees, result) > 0)
        ht <- hash(trees)
        if(!(ht %in% hr)){
          hr <- c(hr, ht)
          result <- c(result, trees)
        }
      }
    }
    if (trace >= 0 &&  (!i%%printevery))
      cat("\rIteration: ", i, ". Best parsimony score so far: ", mp, sep="")
    if ( (kmax >= k) && (i >= minit)) break()
  } # for
  if (trace >= 0)cat("\n")
}  # pratchet
