# Joint reconstruction
joint_pml <- function(x){
  stopifnot(inherits(x, "pml"))
  if(x$k > 1 || x$inv>0) stop("One one rate allowed so far!")
  data <- x$data
  eig <- x$eig
  contrast <- attr(data, "contrast")
  tree <- x$tree
  tree <- reorder(tree, "postorder")
  edge <- tree$edge
  ntip <- Ntip(tree)
  desc <- Descendants(tree, type="children")
  el <- numeric(max(tree$edge))
  el[tree$edge[,2]] <- tree$edge.length
  P <- getP(el * x$rate, eig = eig)
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  levels <- attr(data, "levels")
  allLevels <- attr(data, "allLevels")
  l <- length(tree$edge.length)
  L <- C <- array(NA, c(nr, nc, max(tree$edge)))
  for(i in seq_len(Ntip(tree))){
    L[,,i]<- log(contrast[data[[i]],,drop=FALSE]%*%P[[i]])
  }
  nn <- unique(edge[,1])
  pa <- 0
  root <- edge[l, 1]
  lnr <- seq_len(nr)
  for(i in seq_len(length(nn))){
    pa <- nn[i]
    ch <- desc[[pa]]
    P_i <- P[[pa]]
    tmp1 <- tmp2 <- matrix(0, nr, nc)
    for(j in seq_along(ch)){
      tmp1 <- tmp1 + L[,,ch[j]]
    }
    if(pa==root) break()
    for(j in 1:nc){
      pp <- tmp1 + rep(log(P_i[j,]), each=nr)
      pos <- max.col(pp)
      L[,j,pa]<- pp[ cbind(lnr, pos)]
      C[,j,pa] <- pos
    }
  }
  pp <- tmp1 + rep(log(x$bf), each=nr)
  pos <- max.col(pp)
  L[,,pa] <- pp
  C[,,pa] <- pos
  tree <- reorder(tree)
  if(is.null(tree$node.label)) tree <- makeNodeLabel(tree)
  res <- vector("list", length(tree$node.label))
  names(res) <- tree$node.label
  res[[1]] <- pos
  att <- attributes(data)
  att$names <- tree$node.label
  labels <- c(tree$tip.label, tree$node.label)
  edge <- tree$edge
  nrw <- seq_len(nr)
  for(i in seq_along(edge[,1])){
    ch_i <- edge[i,2]
    pa_i <- edge[i,1]
    if(ch_i > ntip){
      pos <-res[[labels[pa_i]]]
      res[[labels[ch_i]]] <- C[cbind(nrw,pos,ch_i)]
    }
  }
  ind <- match(levels, allLevels)
  for(i in length(res))  res[[i]] <- ind[res[[i]]]
  attributes(res) <- att
  res
}


# Joint reconstruction for parsimony
joint_sankoff <- function(tree, data, cost=NULL){
  stopifnot(inherits(data, "phyDat"))
  stopifnot(inherits(tree, "phylo"))
  tree <- reorder(tree, "postorder")
  edge <- tree$edge
  ntip <- Ntip(tree)
  if (is.null(cost)) {
    levels <- attr(data, "levels")
    l <- length(levels)
    cost <- matrix(1, l, l)
    cost <- cost - diag(l)
  }
  P <- cost

  desc <- Descendants(tree, type="children")
  nr <- attr(data, "nr")
  nc <- attr(data, "nc")
  levels <- attr(data, "levels")
  allLevels <- attr(data, "allLevels")

  l <- nrow(tree$edge)
  C <- array(NA, c(nr, nc, max(tree$edge)))
  tmp <- sankoff(tree, data, cost, "data")
  pscore <- tmp$pscore
  L <- tmp$dat
  nn <- unique(edge[,1])
  pa <- 0
  root <- edge[l, 1]
  lnr <- seq_len(nr)

  for(i in seq_len(length(nn))){
    pa <- nn[i]
    if(pa==root) break()
    for(j in 1:nc){
      pp <- L[[pa]] + rep(P[j,], each=nr)
      pos <- max.col(-pp)
      C[,j,pa] <- pos
    }
  }
  pp <- L[[pa]]
  pos <- max.col(-pp)
  C[,,pa] <- pos
  tree <- reorder(tree)
  if(is.null(tree$node.label)) tree <- makeNodeLabel(tree)
  res <- vector("list", length(tree$node.label))
  names(res) <- tree$node.label
  res[[1]] <- pos
  att <- attributes(data)
  att$names <- tree$node.label
  labels <- c(tree$tip.label, tree$node.label)
  edge <- tree$edge
  nrw <- seq_len(nr)
  for(i in seq_along(edge[,1])){
    ch_i <- edge[i,2]
    pa_i <- edge[i,1]
    if(ch_i > ntip){
      pos <-res[[labels[pa_i]]]
      res[[labels[ch_i]]] <- C[cbind(nrw,pos,ch_i)]
    }
  }
  ind <- match(levels, allLevels)
  for(i in length(res))  res[[i]] <- ind[res[[i]]]
  attributes(res) <- att
  res
}


#' Assign edge length to tree
#'
#' \code{parsimony_edgelength} and \code{acctran} assign edge length to a tree where
#' the edge length is the number of mutations. \code{parsimony_edgelengths}
#' assigns edge lengths using a joint reconstruction based on the sankoff
#' algorithm. Ties are broken at random and trees can be multifurating.
#' \code{acctran} is based on the fitch algorithm and is faster. However trees
#' need to be bifurcating and ties are split.
#' @param tree a tree, i.e. an object of class pml
#' @param data an object of class phyDat
#' @return a tree with edge length.
#' @export
parsimony_edgelength <- function(tree, data){
  if(inherits(tree, "phylo")) return(count_mutations(tree, data))
  if(inherits(tree, "multiPhylo")) {
    res <- lapply(tree, count_mutations, data=data)
    class(res) <- "multiPhylo"
    return(res)
  }
  NULL
}


count_mutations <- function(tree, data){
  tree <- reorder(tree, "postorder")
  data <- data[tree$tip.label]
  tree_tmp <- makeNodeLabel(tree)
  anc <- joint_sankoff(tree_tmp, data)
#  ind <-   length(data)+seq_along(anc)
#  data[ind] <- anc
#  names(dat[ind]) <- names(anc)
  dat <- rbind(data, anc)
  l <- length(dat)
  f <- init_fitch(dat, FALSE, FALSE, m=2L)
  el <- numeric(nrow(tree$edge))
  for(i in seq_along(el)){
    edge_i <-  matrix(c(l+1L, l+1L, tree$edge[i,]), 2, 2)
    el[i] <- f$pscore(edge_i)
  }
  tree$edge.length <- el
  tree
}


map_mutations <- function(tree, data, pos=NULL){
  site <- "sitewise"
  tree <- reorder(tree, "postorder")
  old_tree <- tree
  lev <- attr(data, "allLevels")
  tree <- makeNodeLabel(tree)
  data <- data[c(tree$tip.label, tree$node.label)]
  nr <- attr(data, "nr")
  l <- length(data)
  fun <- function(x, site="pscore", nr){
    if(site=="pscore") return(f$pscore(x))
    sites <- f$sitewise_pscore(x)
    sites[seq_len(nr)]
  }
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  mttns <- vector("list", max(tree$edge))
  index <- attr(data, "index")
  if(is.null(pos)) pos <- seq_along(index)
  M <- as.character(data) # not efficient
  for(i in seq_len(nrow(tree$edge))){
    edge_i <-  matrix(c(l+1L, l+1L, tree$edge[i,]), 2, 2)
    ch_i <- tree$edge[i,2]
    pa_i <- tree$edge[i,1]
    tmp  <- fun(edge_i, site, nr)
    tmp <- which(tmp[index] == 1)
    tmp <- intersect(tmp, pos)
    if(length(tmp > 0)){
      tmp2 <- paste0(M[pa_i, tmp], tmp, M[ch_i, tmp])
      mttns[[ch_i]] <- tmp2
    }
  }
  mttns
}
