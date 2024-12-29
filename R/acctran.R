acctran2 <- function(tree, data) {
  if(!is.binary(tree)) tree <- multi2di(tree)
  tree <- reorder(tree, "postorder")
  edge <- tree$edge
  data <- subset(data, tree$tip.label)
  f <- init_fitch(data, FALSE, FALSE, m=2L)
  psc_node <- f$pscore_node(edge)
  tmp <- reorder(tree)$edge
  tmp <- tmp[tmp[,2]>Ntip(tree), ,drop=FALSE]
  f$traverse(edge)
  if(length(tmp)>0)f$acctran_traverse(tmp)
  psc <- f$pscore_acctran(edge)
  el <- psc
  parent <- unique(edge[,1])
  desc <- Descendants(tree, parent, "children")
  for(i in seq_along(parent)){
    x <- psc_node[parent[i]] -sum(psc[desc[[i]]])
    if(x>0) el[desc[[i]] ] <- el[desc[[i]] ] + x/length(desc[[i]])
  }
  tree$edge.length <- el[edge[,2]]
  tree
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
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pratchet}}
#' @examples
#' data(Laurasiatherian)
#' # lower number of iterations for the example, to run time less than 5 sec.
#' treeRatchet <- pratchet(Laurasiatherian, minit=5, k=5, trace=0)
#' # assign edge length (number of substitutions)
#' treeRatchet <- parsimony_edgelength(treeRatchet, Laurasiatherian)
#' plot(midpoint(treeRatchet))
#' add.scale.bar(0,0, length=100)
#' @export
acctran <- function(tree, data) {
  if (inherits(tree, "multiPhylo")) {
    compress <- FALSE
    if (!is.null(attr(tree, "TipLabel"))){
      compress <- TRUE
      tree <- .uncompressTipLabel(tree)
    }
    res <- lapply(tree, acctran2, data)
    class(res) <- "multiPhylo"
    if (compress) res <- .compressTipLabel(res)
    return(res)
  }
  acctran2(tree, data)
}


#' @rdname acctran
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
