## include in addConfidences, plotBS etc.
transferBootstrap <- function(tree, bstree){
  if(!inherits(bstree, "multiPhylo")) stop("bstrees needs to be of class multiPhylo!")
  bstree <- .uncompressTipLabel(bstree)
  bstree <- .compressTipLabel(bstree, tree$tip.label)
  bstree <- reorder(bstree, "postorder")
  l <- Ntip(tree)
  bp <- prop.part(tree)
  bp <- SHORTwise(bp)[-1]
  not_cherry <- lengths(bp) != 2
  res <- numeric(length(bp))
  for(i in seq_along(bstree)){
     tmp <- bstree[[i]]
     bptmp <- prop.part(tmp)
     bptmp <- SHORTwise(bptmp)[-1]
     ind <- fmatch(bp, bptmp)
     res[!is.na(ind)] <- res[!is.na(ind)] + 1
     # cherries can be check outside
     ind <- which(is.na(ind) & not_cherry)
     for(j in ind) res[j] <- res[j] + Transfer_Index(bp[[j]], tmp$edge, l)
  }
  res <- res / length(bstree) * 100
  tree$node.label <- c(NA_real_, res)
  tree
}

