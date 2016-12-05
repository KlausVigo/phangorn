getAges <- function(x){  
  fun=function(x) max(node.depth.edgelength(x))  
  height=NULL
  if(inherits(x,"phylo")) height <- fun(x)
  if(inherits(x,"multiPhylo")){
    if(!is.null(attr(x, "TipLabel"))){
      x = unclass(x)
      x = .uncompressTipLabel(x)  
      x = unclass(x)  
      height = sapply(x, fun)
    }
    else{
      x = unclass(x)
      height = sapply(x, fun) 
    }
  }
  height
}




#' Plots a densiTree.
#' 
#' An R function to plot trees similar to those produced by DensiTree.
#' 
#' If no consensus tree is provided \code{densiTree} computes a rooted
#' mrp.supertree as a backbone. This should avoid too many unnecessary
#' crossings of edges.  Trees should be rooted, other wise the output may not
#' make sense.
#' 
#' @param x an object of class \code{multiPhylo}.
#' @param type a character string specifying the type of phylogeny, so far
#' "cladogram" (default) or "phylogram" (the default) are supported.
#' @param alpha parameter for semi-transparent colors.
#' @param consensus A tree which is used to define the order of the tip labels.
#' @param optim not yet used.
#' @param scaleX scale trees to have identical heights.
#' @param col edge color.
#' @param width edge width.
#' @param cex a numeric value giving the factor scaling of the tip labels.
#' @param \dots further arguments to be passed to plot.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{plot.phylo}}, \code{\link{plot.networx}}
#' @references densiTree is inspired from the great
#' \href{https://www.cs.auckland.ac.nz/~remco/DensiTree}{DensiTree} program of Remco
#' Bouckaert.
#' 
#' Remco R. Bouckaert (2010) DensiTree: making sense of sets of phylogenetic
#' trees \emph{Bioinformatics}, \bold{26 (10)}, 1372-1373.
#' @keywords plot
#' @examples
#'   
#' data(Laurasiatherian)
#' set.seed(1)
#' bs <- bootstrap.phyDat(Laurasiatherian, FUN = 
#'    function(x)upgma(dist.hamming(x)), bs=25)
#' # cladogram nice to show topological differences
#' densiTree(bs, optim=TRUE, type="cladogram", col="blue")
#' densiTree(bs, optim=TRUE, type="phylogram", col="green")
#' \dontrun{
#' # phylograms are nice to show different age estimates
#' require(PhyloOrchard)
#' data(BinindaEmondsEtAl2007)
#' BinindaEmondsEtAl2007 <- .compressTipLabel(BinindaEmondsEtAl2007) 
#' densiTree(BinindaEmondsEtAl2007, type="phylogram", col="red")
#' }
#' 
#' 
#' @export densiTree
densiTree <- function(x, type="cladogram", alpha=1/length(x), consensus=NULL, optim=FALSE, scaleX=FALSE, col=1, width=1, cex=.8, ...) {
  if(!inherits(x,"multiPhylo"))stop("x must be of class multiPhylo")
  compressed <- ifelse(is.null(attr(x, "TipLabel")), FALSE, TRUE)
  if(is.null(consensus))consensus <- superTree(x)
  if(inherits(consensus,"multiPhylo")) consensus = consensus[[1]]
  type <- match.arg(type, c("phylogram", "cladogram"))
  consensus = reorder(consensus, "postorder")
  e2 = reorder(consensus)$edge[,2]
  nTip = as.integer(length(consensus$tip.label))
  tiporder = e2[e2<=nTip]   
  maxBT = max(getAges(x))
  if(scaleX) maxBT=1.0
  label = rev(pretty(c(maxBT,0)))
  maxBT = max(label)
  xy = plotPhyloCoor(consensus, ...)
  yy = xy[,2]
  plot.new() 
  tl = which.max(nchar(consensus$tip.label))
  sw <- strwidth(consensus$tip.label[tl],cex=cex) * 1.1
  plot.window(xlim=c(0, 1.0+sw), ylim=c(0, nTip+1))
  axis(side=1,at=seq(0,1.0, length.out=length(label)), labels=label)
  text(x=rep(1.0,Ntip(consensus)),y=yy[1:nTip],labels=consensus$tip.label,pos=4,cex=cex)  
  tip.order = yy[1:nTip]
  for (treeindex in 1:length(x)) {
    tmp <- reorder(x[[treeindex]], "postorder")
    if(!compressed) tip.order <- match(tmp$tip.label, consensus$tip.label)
    xy <- plotPhyloCoor(tmp, tip.order=tiporder, ...)
    xx = xy[,1]
    yy = xy[,2]
    if(scaleX) xx <- xx/max(xx)
    else xx <- xx/maxBT 
    xx <- xx + (1.0 - max(xx))
    e1=tmp$edge[,1]
    e2=tmp$edge[,2]
    if(type=="cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color=adjustcolor(col, alpha.f=alpha), edge.width=width, edge.lty=1)
    if(type=="phylogram"){
      Ntip <- min(e1)-1L 
      Nnode <- tmp$Nnode 
      phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, TRUE, edge.color=adjustcolor(col, alpha.f=alpha), edge.width=width, 1) 
    }
  }  
}






