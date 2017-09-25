nodeHeight <- function(tree) 
{
    if(is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise")
        tree <- reorder(tree, "postorder")
    edge <- tree$edge[, 2]
    node <- tree$edge[, 1]
    m <- max(tree$edge)
    el <- double(m)
    el[edge] <- tree$edge.length
    res <- .C("nodeH", as.integer(edge), as.integer(node), el, 
              as.integer(length(edge)), double(m))[[5]] 
    max(res) - res
}



ancstat <- function(phy, x){
  contrast <- attr(x, "contrast")
  storage.mode(contrast) <- "integer"
  phy <- reorder(phy, "postorder")
  res <- matrix(0L, max(phy$edge), ncol(contrast))
  colnames(res) <- attr(x,"levels")
  nTips <- length(phy$tip.label)
  pa <- phy$edge[,1]
  ch <- phy$edge[,2]
  res[1:nTips, ] <- contrast[as.numeric(x)[match(phy$tip.label, names(x))],, 
            drop=FALSE]
  for(i in seq_along(pa)){
    res[pa[i],] <- res[pa[i],] | res[ch[i],]    
  }
  res
}


comp <- function(x, y){
  tmp1 <- matrix(rowSums(x), nrow(x), nrow(y))
  res <- matrix(rowSums(y), nrow(x), nrow(y), byrow=TRUE)
  tmp3 <- tcrossprod(x, 1-y)  
  tmp0 <- tcrossprod(x, y)
  tmp0[tmp3>0] <- 0L
  res[!(tmp0>(tmp1 - 1e-8))] <- 10000000L 
  apply(res, 1, which.min)
}


comp2 <- function(x, y){
  res <- matrix(rowSums(x), nrow(x), nrow(y))
  tmp1 <- matrix(rowSums(y), nrow(x), nrow(y), byrow=TRUE)
  tmp3 <- tcrossprod(1-x, y)  
  tmp0 <- tcrossprod(x, y)
  tmp0[tmp3>0] <- 0L
  res[tmp0<2] <- Inf
  apply(res, 2, which.min)
}

# single linkage of minimal coalescent times
# extends speciesTree fom ape



#' Species Tree
#' 
#' \code{coalSpeciesTree} estimates species trees and can handle multiple
#' individuals per species.
#' 
#' \code{coalSpeciesTree} estimates a single linkage tree as suggested by Liu
#' et al. (2010) from the element wise minima of the cophenetic matrices of the
#' gene trees. It extends \code{speciesTree} in ape as it allows that have
#' several individuals per gene tree.
#' 
#' @export
#' 
#' @param tree an object of class \code{multiPhylo}
#' @param X A \code{phyDat} object to define which individual belongs to which
#' species.
#' @param sTree A species tree which fixes the topology.
#' @return The function returns an object of class \code{phylo}.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} Emmanuel Paradies
#' @seealso \code{\link{speciesTree}}
#' @references Liu, L., Yu, L. and Pearl, D. K. (2010) Maximum tree: a
#' consistent estimator of the species tree. \emph{Journal of Mathematical
#' Biology}, \bold{60}, 95--106.
#' @keywords cluster
#' 
coalSpeciesTree <- function(tree, X=NULL, sTree=NULL){
  
  if(is.null(X))return(speciesTree(tree))  
  trees <- unclass(tree)
  States <- lapply(tree, ancstat, X)
  NH <- lapply(tree, nodeHeight)
  if(is.null(sTree)){
    l <- attr(X, "nc")
    m <- choose(l, 2)
    SST <- matrix(0L, m, l)
    k <- 1
    for(i in 1:(l-1)){
      for(j in (i+1):l){
        SST[k, i] <- SST[k,j] <- 1L
        k <- k+1
      }
    }
    Y <- matrix(Inf, length(NH), nrow(SST)) 
    dm <- rep(Inf, m)
    for(i in seq_along(NH)){
      ind <- comp2(States[[i]],SST)  
      dm <- pmin(dm, NH[[i]][ind])
      #   for(j in 1:length(ind))Y[i, ind[j]] = min(Y[i, ind[j]], NH[[i]][j])
    }
    dm <- structure(2*dm, Labels = attr(X, "levels"), Size = l, class = "dist", 
                    Diag = FALSE, Upper = FALSE)
    
    sTree <- upgma(dm, "single")   
    # dm of pairwise states
  }   
  else{ 
    SST <- ancstat(sTree, X)
    Y <- matrix(Inf, length(NH), nrow(SST)) 
    for(i in seq_along(NH)){
      ind <- comp(States[[i]],SST) 
      for(j in seq_along(ind))Y[i, ind[j]] <- min(Y[i, ind[j]], NH[[i]][j])
    }
    STH <- apply(Y, 2, min)
    sTree$edge.length <- STH[sTree$edge[,1]] - STH[sTree$edge[,2]]
  }
  sTree
}

