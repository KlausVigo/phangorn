#
# Maximum Parsimony
#
sankoff.quartet <- function(dat, cost, p, l, weight) {
  erg <- .Call('sankoffQuartet', sdat = dat, sn = p, scost = cost, sk = l)
  sum(weight * erg)
}


#' Parsimony tree.
#'
#'
#' \code{parsimony} returns the parsimony score of a tree using either the
#' sankoff or the fitch algorithm. \code{optim.parsimony} tries to find the
#' maximum parsimony tree using either Nearest Neighbor Interchange (NNI)
#' rearrangements or sub tree pruning and regrafting (SPR). \code{pratchet}
#' implements the parsimony ratchet (Nixon, 1999) and is the preferred way to
#' search for the best tree.  \code{random.addition} can be used to produce
#' starting trees.
#'
#' The "SPR" rearrangements are so far only available for the "fitch" method,
#' "sankoff" only uses "NNI". The "fitch" algorithm only works correct for
#' binary trees.
#'
#' @aliases parsimony
#' @aliases optim.parsimony sankoff fitch pratchet
#' random.addition acctran
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
  attrData$weight <- tapply(attrData$weight, index, sum)
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
#' @param tree tree to start the nni search from.
#' @param cost A cost matrix for the transitions between two states.
#' @param sitewise return CI/RI for alignment or sitewise
#'
#' @seealso \code{\link{parsimony}}, \code{\link{pratchet}},
#' \code{\link{fitch}}, \code{\link{sankoff}}, \code{\link{bab}},
#' \code{\link{ancestral.pars}}
#'
#' @rdname CI
#' @export
CI <- function(tree, data, cost = NULL, sitewise = FALSE) {
  if (sitewise) pscore <- sankoff(tree, data, cost = cost, site = "site")
  else pscore <- sankoff(tree, data, cost = cost)
  weight <- attr(data, "weight")
  data <- subset(data, tree$tip.label)
  m <- lowerBound(data, cost = cost)
  if (sitewise) {
    return( (m / pscore)[attr(data, "index")])
  }
  sum(m * weight) / pscore
}


#' @rdname CI
#' @export
RI <- function(tree, data, cost = NULL, sitewise = FALSE) {
  if (sitewise) pscore <- sankoff(tree, data, cost = cost, site = "site")
  else pscore <- sankoff(tree, data, cost = cost)
  data <- subset(data, tree$tip.label)
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



#
# Sankoff
#


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


prepareDataSankoff <- function(data) {
  contrast <- attr(data, "contrast")
  contrast[contrast == 0] <- 1e+06
  contrast[contrast == 1] <- 0
  for (i in seq_along(data)) data[[i]] <- contrast[data[[i]], , drop = FALSE]
  data
}


fit.sankoff <- function(tree, data, cost,
                        returnData = c("pscore", "site", "data")) {
  tree <- reorder(tree, "postorder")
  returnData <- match.arg(returnData)
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  weight <- attr(data, "weight")
  nr <- attr(data, "nr")
  q <- length(tree$tip.label)
  nc <- attr(data, "nc")
  m <- length(edge) + 1
  dat <- vector(mode = "list", length = m)
  dat[1:q] <- data[tree$tip.label]
  node <- as.integer(node - 1)
  edge <- as.integer(edge - 1)
  mNodes <- as.integer(max(node) + 1)
  tips <- as.integer( (seq_along(tree$tip.label)) - 1)
  res <- .Call('sankoff3', dat, as.numeric(cost), as.integer(nr),
    as.integer(nc), node, edge, mNodes, tips)
  root <- getRoot(tree)
  erg <- .Call('C_rowMin', res[[root]], as.integer(nr), as.integer(nc))
  if (returnData == "site") return(erg)
  pscore <- sum(weight * erg)
  result <- pscore
  if (returnData == "data") {
    result <- list(pscore = pscore, dat = res)
  }
  result
}


pnodes <- function(tree, data, cost) {
  tree <- reorder(tree, "postorder")
  node <- tree$edge[, 1]
  edge <- tree$edge[, 2]
  nr <- nrow(data[[1]])
  nc <- ncol(data[[1]])
  node <- as.integer(node - 1)
  edge <- as.integer(edge - 1)
  .Call('pNodes', data, as.numeric(cost), as.integer(nr), as.integer(nc),
    node, edge)
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
  # wahrscheinlich schneller: cvector <- allCildren(tree)
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


sankoff.nni <- function(tree, data, cost, ...) {
  if (is.rooted(tree)) tree <- reorder(unroot(tree), "postorder")
  INDEX <-  indexNNI(tree)
  rootEdges <- attr(INDEX, "root")
  if (!inherits(data, "phyDat"))
    stop("data must be of class phyDat")
  levels <- attr(data, "levels")
  l <- length(levels)
  weight <- attr(data, "weight")
  p <- attr(data, "nr")
  i <- 1
  tmp <- fit.sankoff(tree, data, cost, returnData = "data")
  p0 <- tmp[[1]]
  datf <- tmp[[2]]
  datp <- pnodes(tree, datf, cost)

  parent <- tree$edge[, 1]
  m <- dim(INDEX)[1]
  k <- min(parent)
  pscore <- numeric(2 * m)

  for (i in 1:m) {
    ei <- INDEX[i, ]
    datn <- datf[ei[1:4]]
    if (!(ei[5] %in% rootEdges)) datn[1] <- datp[ei[1]]
    pscore[(2 * i) - 1] <- sankoff.quartet(datn[ c(1, 3, 2, 4)],
      cost, p, l, weight)
    pscore[(2 * i)] <- sankoff.quartet(datn[ c(1, 4, 3, 2)],
      cost, p, l, weight)
  }
  swap <- 0
  candidates <- pscore < p0
  while (any(candidates)) {

    ind <- which.min(pscore)
    pscore[ind] <- Inf
    if (ind %% 2) swap.edge <- c(2, 3)
    else swap.edge <- c(2, 4)

    tree2 <- changeEdge(tree, INDEX[(ind + 1) %/% 2, swap.edge])
    test <- fit.sankoff(tree2, data, cost, "pscore")

    if (test >= p0) candidates[ind] <- FALSE
    if (test < p0) {
      p0 <- test
      swap <- swap + 1
      tree <- tree2
      candidates[ind] <- FALSE
      indi <- which(rep(colSums(apply(INDEX, 1, match, INDEX[(ind + 1) %/% 2, ],
        nomatch = 0)) > 0, each = 2))
      candidates[indi] <- FALSE
      pscore[indi] <- Inf
    }
  }
  list(tree = tree, pscore = p0, swap = swap)
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
  eps <- 1e-08
  trace <- trace - 1

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
    return(tree)
  }

  if (perturbation != "random_addition"){
    if(is.null(start)) start <- optim.parsimony(fastme.ols(dist.hamming(data)),
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
      print(paste("Best pscore so far:", attr(tree, "pscore")))
  }
  FUN <- function(data, tree, method, rearrangements, ...)
    optim.parsimony(tree, data = data, method = method,
                    rearrangements = rearrangements, ...)
  result <- tree
  if(!is.null(attr(data, "duplicated"))){
    result <- addTaxa(result, attr(data, "duplicated"))
  }
  on.exit({
    if (!all && inherits(result, "multiPhylo")) result <- result[[1]]
#    if(!is.null(attr(data, "duplicated")))
#      result <- addTaxa(result, attr(data, "duplicated"))
    #    else class(result) <- "multiPhylo"
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
      if(length(bs_ind) > 0)p_trees <- optim.parsimony(tree, bs_data,
          trace = trace, method = method, rearrangements = rearrangements, ...)
      else p_trees <- stree(length(data), tip.label = names(data))
      trees <- optim.parsimony(p_trees, data, trace = trace,
                     method = method, rearrangements = rearrangements, ...)
    }
    if (perturbation == "stochastic") {
      p_trees <- rNNI(tree, floor(nTips / 2))
      trees <- optim.parsimony(p_trees, data, trace = trace, method = method,
                               rearrangements = rearrangements, ...)
    }
    if (perturbation == "random_addition") {
      p_trees <- random.addition(data)
      trees <- optim.parsimony(p_trees, data, trace = trace, method = method,
                               rearrangements = rearrangements, ...)
    }
    if(!is.null(attr(data, "duplicated"))){
      p_trees <- addTaxa(p_trees, attr(data, "duplicated"))
      trees <- addTaxa(trees, attr(data, "duplicated"))
    }
    start_trees[[i]] <- p_trees
    search_trees[[i]] <- trees
    pscores <- attr(trees, "pscore")
    mp1 <- min(pscores)
    if ( (mp1 + eps) < mp) {
      kmax <- 1
      result <- trees
      tree <- trees
      mp <- mp1
    }
    else{
      kmax <- kmax + 1
      if( all && (mp1 < (mp + eps)) && all(RF.dist(trees, result) > 0))
        result <- c(result, trees)
    }
    if (trace >= 0)
      print(paste("Best pscore so far:", mp))
    if ( (kmax >= k) && (i >= minit)) break()
  } # for
}  # pratchet


optim.sankoff <- function(tree, data, cost = NULL, trace = 1, ...) {
  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  if (is.rooted(tree)) tree <- unroot(tree)
  tree <- reorder(tree, "postorder")
  if (!inherits(data, "phyDat")) stop("data must be of class phyDat")
  addTaxa <- FALSE
  mapping <- map_duplicates(data)
  if (!is.null(mapping)) {
    addTaxa <- TRUE
    tree2 <- drop.tip(tree, mapping[, 1])
    tree2 <- unroot(tree2)
    tree <- reorder(tree2, "postorder")
  }

  rt <- FALSE
  dat <- prepareDataSankoff(data)
  l <- attr(dat, "nc")
  if (is.null(cost)) {
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }


  tree$edge.length <- NULL
  swap <- 0
  iter <- TRUE
  pscore <- fit.sankoff(tree, dat, cost, "pscore")

  on.exit({
    if (rt) tree <- acctran(tree, data)
    if (addTaxa) {
      if (rt) tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2],
          edge.length = rep(0, nrow(mapping)))
      else tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2])
    }
    attr(tree, "pscore") <- pscore
    return(tree)
  })

  while (iter) {
    res <- sankoff.nni(tree, dat, cost, ...)
    tree <- res$tree
    if (trace > 1) cat("optimize topology: ", pscore, "-->", res$pscore, "\n")
    pscore <- res$pscore
    swap <- swap + res$swap
    if (res$swap == 0) iter <- FALSE
  }
  if (trace > 0) cat("Final p-score", pscore, "after ", swap,
                     "nni operations \n")
}
