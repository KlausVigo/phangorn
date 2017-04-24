#' Maximum agreement subtree
#' 
#' \code{mast} computes the maximum agreement subtree (MAST). 
#' 
#' The code is derived from exames in Valverde (2009), for the original code see 
#' \url(http://www.cs.upc.edu/~valiente/comput-biol/). 
#' 
#' @param x a tree, i.e. an object of class \code{phylo}.
#' @param y a tree, i.e. an object of class \code{phylo}.
#' @param tree a logicical ,if TRUE returns a tree other wise the tip labels 
#' of the the maximum agreement subtree.
#' @return \code{mast} returns a vector of the tip labels in the MAST.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com} 
#' @seealso \code{\link{SPR.dist}}
#' @references
#' G. Valiente (2009). \emph{Combinatorial Pattern Matching Algorithms in Computational Biology using Perl and R}. Taylor & Francis/CRC Press 
#' 
#'
#' @keywords cluster
#' @examples
#' tree1 <- rtree(100, rooted=FALSE)
#' tree2 <- rSPR(tree1, 5)
#' tips <- mast(tree1, tree2)
#' 
#' @rdname mast
#' @export
mast <- function (x, y, tree=TRUE) {
    # make x and y have same taxa set
    shared_tips <- intersect(x$tip.label, y$tip.label)
    if(length(shared_tips) < length(x$tip.label)) 
        x <- drop.tip(x, setdiff(x$tip.label, shared_tips))
    if(length(shared_tips) < length(y$tip.label)) 
        y <- drop.tip(y, setdiff(y$tip.label, shared_tips))
    
    is.ancestor <- function(anc, des, pvec) {
        if (anc==des) return(TRUE)
        return(anc %in% pvec)
    }
    
    # make order of labels the same
    y <- .compressTipLabel(c(y), x$tip.label)[[1]]
    x <- reorder(x, "postorder")
    y <- reorder(y, "postorder")
    
#    po1 <- postorder(x)
#    po2 <- postorder(y)
 
    po1 <- c(x$edge[,2], x$edge[nrow(x$edge),1])
    po2 <- c(y$edge[,2], y$edge[nrow(y$edge),1])
    
       
    nTips <- length(x$tip.label)
    
    p_vec_1 <- Ancestors(x, 1L:max(x$edge))
    p_vec_2 <- Ancestors(y, 1L:max(y$edge))
    
    CH1 <- allChildren(x)
    CH2 <- allChildren(y)
    
    m <- matrix(rep(list(),length(po1)),nrow=length(po1),ncol=length(po2),
                dimnames=list(po1,po2))
    for (i in po1) {
        for (j in po2) {
            mm <- list()
            if (i<=nTips) {
                if (j<=nTips) {
                    if (i==j)
                        mm <- c(x$tip.label[i])
                } else {
                    if(is.ancestor(j, i, p_vec_2[[i]]))
                        mm <- c(x$tip.label[i])
                }
            } else {
                if (j <= nTips) {
                    if(is.ancestor(i, j, p_vec_1[[j]]))
                        mm <- c(y$tip.label[j])
                } else {
                    l1 <- CH1[[i]][1]
                    r1 <- CH1[[i]][2]
                    l2 <- CH2[[j]][1]
                    r2 <- CH2[[j]][2]
                    mm <- c(m[[l1,l2]],m[[r1,r2]])
                    if ( length(m[[l1,r2]]) + length(m[[r1,l2]]) > length(mm) )
                        mm <- c(m[[l1,r2]],m[[r1,l2]])
                    if ( length(m[[i,l2]]) > length(mm) )
                        mm <- m[[i,l2]]
                    if ( length(m[[i,r2]]) > length(mm) )
                        mm <- m[[i,r2]]
                    if ( length(m[[l1,j]]) > length(mm) )
                        mm <- m[[l1,j]]
                    if ( length(m[[r1,j]]) > length(mm) )
                        mm <- m[[r1,j]]
                }
            }
            m[[i,j]] <- mm
        }
    }
    res <- unlist(m[[length(x$tip.label)+1L,length(y$tip.label)+1L]])
    if(tree) res <- drop.tip(x, setdiff(x$tip.label, res))
    res
}
