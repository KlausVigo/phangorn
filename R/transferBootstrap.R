## include in addConfidences, plotBS etc.
transferBootstrap <- function(tree, bstree){
  if(!inherits(bstree, "multiPhylo")) stop("bstrees needs to be of class multiPhylo!")
  bstree <- .uncompressTipLabel(bstree)
  bstree <- .compressTipLabel(bstree, tree$tip.label)
  bstree <- reorder(bstree, "postorder")
  l <- Ntip(tree)
  bp <- bipart(tree)
  bp <- SHORTwise(bp, l, TRUE)
  not_cherry <- lengths(bp) != 2
  res <- numeric(length(bp))
  for(i in seq_along(bstree)){
     tmp <- bstree[[i]]
     bptmp <- bipart(tmp)
     bptmp <- SHORTwise(bptmp, l, TRUE)
     ind <- fmatch(bp, bptmp)
     res[!is.na(ind)] <- res[!is.na(ind)] + 1
     # cherries can be check outside
     ind <- which(is.na(ind) & not_cherry)
     for(i in ind) res[i] <- res[i] + Transfer_Index(bp[[i]], tmp$edge, l)
  }
  res <- res / length(bstree)
  tree$node.label <- c(NA_real_, res)
  tree
}


# R version (C++ ~ 10x faster)
transfer_index <- function(bp, edge, l){
#  edge <- tree$edge
  m <- max(edge)
  p <- length(bp)
  lmp <- l - p
  best <- p - 1
  l0 <- l1 <- numeric(m)
  l1[bp] <- 1
  l0[seq_len(l)[-bp]] <- 1
  node <- edge[1,1]
  for(i in seq_len(nrow(edge))){
    ni <- edge[i, 1]
    ei <- edge[i, 2]
    l0[ni] <- l0[ni] + l0[ei]
    l1[ni] <- l1[ni] + l1[ei]
    if(ni != node){
      best <- min(best, (p - l1[node]) + l0[node], (lmp - l0[node]) + l1[node])
      if(best==1) return(1 - (best / (p-1)))
      node <- ni
    }
  }
  tmp <- min((p - l1[node]) + l0[node], (lmp - l0[node]) + l1[node])
  best <- min(best, tmp)
  return(1 - (best / (p-1)))
}

