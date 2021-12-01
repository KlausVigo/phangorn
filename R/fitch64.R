init_fitch <- function(obj, parsinfo=FALSE, order=FALSE, m=4L, ...){
  if(parsinfo) obj <- removeParsimonyUninfomativeSites(obj, ...)
  if(is.null(attr(obj, "p0"))) attr(obj, "p0") <- 0
  attr(obj, "nSeq") <- length(obj) # add lengths
  if(is.null(attr(obj, "weight"))){
    attr(obj, "weight") <- rep(1, length(obj[[1]]))
    order <- FALSE
  }
  if(sum( abs (attr(obj, "weight") %% 1L) ) >1) order <- FALSE
  if(order){
    ord <- order(attr(obj, "weight"), decreasing = TRUE)
    obj <- subset(obj, select = ord) # inside C (1 less copy)
  }
  weight <- attr(obj, "weight")
  l <- length(weight)
  # can be NA
  first_1 <- match(1L, weight)
  if(is.na(first_1)) first_1 <- l
  if(first_1 == 1L) first_1 <- 0L
  if(!all(weight[first_1:l]==1)) first_1 <- l

  d_con <- dim(attr(obj, "contrast"))
  contrast <- matrix(0L, d_con[1], d_con[2])
  contrast[attr(obj, "contrast") > 1e-8] <- 1L
  contrast <- rbind(0L, contrast)
  storage.mode(contrast) <- "integer"
  attr(obj, "contrast") <- contrast
  f <- new(Fitch, obj, as.integer(first_1), as.integer(m))
  f
}


#' @rdname parsimony
#' @export
fitch <- function(tree, data, site = "pscore"){
  if(any(!is.binary(tree))) tree <- multi2di(tree)
  tree <- reorder(tree, "postorder")
  nr <- attr(data, "nr")
  fun <- function(tree, site="pscore", nr){
    if(site=="pscore") return(f$pscore(tree$edge))
    sites <- f$sitewise_pscore(tree$edge)
    sites[seq_len(nr)]
  }
  fun2 <- function(tree, data, site, nr) {
    data <- subset(data, tree$tip.label)
    f <- init_fitch(data, FALSE, FALSE, m=2L)
    if(site=="pscore") return(f$pscore(tree$edge))
    sites <- f$sitewise_pscore(tree$edge)
    sites[seq_len(nr)]
  }
  if (inherits(tree, "multiPhylo")) {
    TL <- attr(tree, "TipLabel")
    if (!is.null(TL)) {
      data <- subset(data, TL)
      f <- init_fitch(data, FALSE, FALSE, m=2L)
      tree <- unclass(tree)
      res <- sapply(tree, function(x)fun(x, site=site, nr=nr))
    }
    else{
      res <- sapply(tree, fun2, data, site, nr)
    }
    return(res)
  }
  if(inherits(tree, "phylo")) {
    data <- subset(data, tree$tip.label)
    f <- init_fitch(data, FALSE, FALSE, m=2L)
    return(fun(tree, site, nr))
  }
  NULL
}



#' @rdname parsimony
#' @export
random.addition <- function (data, tree=NULL, method = "fitch")
{
  label <- names(data)
  nTips <- as.integer(length(label))
  if (nTips < 4L)
    return(stree(nTips, tip.label = sample(label)))
  if(!is.null(tree)){
    if(!is.binary(tree)) tree <- multi2di(tree)
    if(!is.null(tree$edge.length)) tree$edge.length <- NULL
    tree <- reorder(tree, "postorder")
    tips <- Ntip(tree)
    edge <- tree$edge
    edge[edge>tips] <- edge[edge>tips] + as.integer(nTips - tips)
    tree$edge <- edge
    remaining <- sample(setdiff(label, tree$tip.label))
    tree$tip.label <- c(tree$tip.label, remaining)
    tree <- checkLabels(tree, label)
    remaining <- match(remaining, label)
  }

  else{
    remaining <- as.integer(sample(nTips))
    tree <- structure(list(edge = structure(c(rep(nTips + 1L, 3),
                        remaining[1:3]), .Dim = c(3L, 2L)), tip.label = label,
                        Nnode = 1L), .Names = c("edge", "tip.label", "Nnode"),
                        class = "phylo", order = "postorder")
    remaining <- remaining[-c(1:3)]
  }
  f <- init_fitch(data, parsinfo = TRUE, order = TRUE, m=4L)
  for (i in remaining) {
    edge <- tree$edge
    f$traversetwice(edge, 0L)
    f$root_all_node(edge)
    score <- f$pscore_vec(edge[,2] + 2 * nTips, i)
    nt <- which.min(score)
    tree <- addOne(tree, i, nt)
  }
  attr(tree, "pscore") <- f$pscore(tree$edge)
  tree
}


fitch_spr <- function (tree, f, trace=0L)
{
  nTips <- as.integer(length(tree$tip.label))
  m <- max(tree$edge)
#  f <- init_fitch(data, FALSE, FALSE, m=4L)
  for (i in 1:nTips) {
# remove tip
    treetmp <- dropTip(tree, i)
    edge <- treetmp$edge
    f$prep_spr(edge)
    score <- f$pscore_vec(edge[,2] + 2 * nTips, i)
    nt <- which.min(score)
# check if different
    tree <- addOne(treetmp, i, nt)
  }
  root <- getRoot(tree)
  ch <- allChildren(tree)
  for (i in (nTips + 1L):m) {
    if (i != root) {
      tmp <- dropNode(tree, i, all.ch = ch)
      if (!is.null(tmp)) {
        f$prep_spr(tmp[[1]]$edge)
        score <- f$pscore_vec(tmp[[1]]$edge[,2] + 2 * nTips, i)
        nt <- which.min(score)
        if(!(tmp[[1]]$edge[nt, 2L] %in% tmp[[4]])){
          tree <- addOneTree(tmp[[1]], tmp[[2]], nt, tmp[[3]])
          ch <- allChildren(tree)
          if(trace) print(f$pscore(tree$edge))
        }
      }
    }
  }
  attr(tree, "pscore") <- f$pscore(tree$edge)
  tree
}


indexNNI_fitch <- function(tree) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  ind <- child
  nTips <- length(tree$tip.label)
  ind <- ind[ind > nTips]
  edgeMatrix <- matrix(0L, length(ind), 6L)

  pvector <- integer(max(parent))
  pvector[child] <- parent
  cvector <- Children(tree) # allChildren
  #     a         d
  #      \       /
  #       e-----f       d is closest to root, f is root from subtree a,b,c
  #      /       \
  #     b         c     c(a,b,c,d,e,f)

  k <- 1
  for (i in ind) {
    f <- pvector[i]
    ab <- cvector[[i]]
    ind1 <- cvector[[f]]
    cd <- ind1[ind1 != i]
    if (pvector[f]){
      cd <- c(cd, f + 2L * nTips)
      ef <- c(i, f)
    }
    else ef <- c(i, cd[2])
    edgeMatrix[k, ] <- c(ab, cd, ef)
    k <- k + 1
  }
  #cbind(edgeMatrix[c(1, 3, 2, 4), ], edgeMatrix[c(2, 3, 1, 4), ])
  edgeMatrix
}


nni2 <- function(x){
  INDEX <- indexNNI_fitch(x)[, 1:4]
  INDEX <- rbind(INDEX[, c(1, 3, 2, 4)], INDEX[, c(2, 3, 1, 4)])
  l <- nrow(INDEX)
  res <- vector("list", l)
  #  for(i in seq_len(l)) res[[i]] <- changeEdge(x, INDEX[c(2, 3), i])
  for(i in seq_len(l)) res[[i]] <- changeEdge(x, INDEX[i, c(2, 3)])
  class(res) <- "multiPhylo"
  res
}


fitch_nni <- function(tree, f) {
  p0 <- f$pscore(tree$edge)
  nTips <- as.integer(length(tree$tip.label))
  INDEX <- indexNNI_fitch(tree)
  l <- nrow(INDEX)
  f$traversetwice(tree$edge, 1L)
  M <- f$pscore_nni(INDEX[, 1L:4L, drop=FALSE])
  M <- M[, -1L] - M[, 1L]
  M <- as.vector(M)
  INDEX <- rbind(INDEX[, c(1, 3, 2, 4, 5, 6)], INDEX[, c(2, 3, 1, 4, 5, 6)])
  swap <- 0
  candidates <- which(M < 0)
  while (length(candidates)>0) {
    pscore <- M[candidates]
    ind <- which.min(pscore)
    tree2 <- changeEdge(tree, INDEX[candidates[ind], c(2, 3)])
    test <- f$pscore(tree2$edge)
    if (test >= p0)
      candidates <- candidates[-ind]
    if (test < p0) {
      p0 <- test
      swap <- swap + 1
      tree <- tree2
      indi <- which(INDEX[, 5] %in% INDEX[candidates[ind],])
      candidates <- setdiff(candidates, indi)
    }
  }
  p0 <- f$pscore(tree$edge)
  list(tree = tree, pscore = p0, swap = swap)
}


optim.fitch <- function(tree, data, trace = 1, rearrangements = "NNI", ...) {
  if (!inherits(tree, "phylo")) stop("tree must be of class phylo")
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
    attr(tree, "order") <- NULL
  }
  if (is.rooted(tree)) {
    tree <- unroot(tree)
    attr(tree, "order") <- NULL
  }
  if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
    tree <- reorder(tree, "postorder")
  if (class(data)[1] != "phyDat") stop("data must be of class phyDat")

  rt <- FALSE

  #  New
  data <- removeParsimonyUninfomativeSites(data, recursive=TRUE)

  star_tree <- ifelse(attr(data, "nr") == 0, TRUE, FALSE)
  add_taxa <- ifelse(is.null(attr(data, "duplicated")), FALSE, TRUE)
  nTips <- length(data)
  if (nTips < 4L  || star_tree) {
    nam <- names(data)
    if (star_tree) tree <- stree(length(nam), tip.label = nam)
    else tree <- stree(nTips, tip.label = nam)
    if(add_taxa) tree <- addTaxa(tree, attr(data, "duplicated"))
    tree <- unroot(tree)
    return(tree)
  }

  tree <- keep.tip(tree, names(data))
  if(length(tree$tip.label) > 2) tree <- unroot(tree)
  tree <- reorder(tree, "postorder")
  p0 <- attr(data, "p0")

  nr <- attr(data, "nr")
  nTips <- as.integer(length(tree$tip.label))
  if(nTips < 5) rearrangements <- "NNI"

  data <- subset(data, tree$tip.label, order(attr(data, "weight"),
                                             decreasing = TRUE))
  f <- init_fitch(data, FALSE, FALSE, m=4L)

  m <- nr * (2L * nTips - 2L)
  on.exit({
    if (add_taxa) tree <- addTaxa(tree, attr(data, "duplicated"))
    tree <- unroot(tree)
    attr(tree, "pscore") <- pscore
    return(tree)
  })

  tree$edge.length <- NULL
  swap <- 0
  iter <- TRUE
  if(nTips < 4) iter <- FALSE
  pscore <- f$pscore(tree$edge)
  while (iter) {
    res <- fitch_nni(tree, f)
    tree <- res$tree
    psc <- res$pscore
    if (trace > 1) cat("optimize topology (NNI): ", pscore, "-->", psc, "\n")
    if(psc < pscore) pscore <- psc
    swap <- swap + res$swap
    if (res$swap == 0) {
      if (rearrangements == "SPR") {
        tree2 <- fitch_spr(tree, f)
        psc <- f$pscore(tree2$edge)
        if (trace > 1) cat("optimize topology (SPR): ", pscore, "-->",
                           psc , "\n")
        if (pscore < psc + 1e-6) iter <- FALSE
        else{
          pscore <- psc
          tree <- tree2
        }
      }
      #      if (rearrangements == "TBR") {}
      else iter <- FALSE
    }
  }
  if (trace > 0) cat("Final p-score", pscore, "after ", swap,
                     "nni operations \n")
}

