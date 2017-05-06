getAges <- function(x){  
  fun=function(x) max(node.depth.edgelength(x))  
  height=NULL
  if(inherits(x,"phylo")) height <- fun(x)
  if(inherits(x,"multiPhylo")){
    if(!is.null(attr(x, "TipLabel"))){
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


add_tiplabels <- function(xy, tip.label, direction, adj, font, srt=0, cex=1, 
        col=1, label.offset=0){
    direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards", 
        "downwards"))
    horizontal <- direction %in% c("rightwards", "leftwards")
    nTips <- length(tip.label)
    xx <- rep(1, nrow(xy))
    yy <- xy[,2 ]
    if(direction == "leftwards" | direction == "downwards") xx = xx*0
    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- xy[,1]
    }
    MAXSTRING <- max(strwidth(tip.label, cex = cex))
    loy <- 0
    if(direction == "rightwards") lox <- label.offset + MAXSTRING * 1.05 * adj
    if(direction == "leftwards") lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
    if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
            loy <- -loy
            srt <- 180 + srt
        }
    }
    text(xx[1:nTips] + lox, yy[1:nTips] + loy, tip.label, adj = adj, 
         font = font, srt = srt, cex = cex, col = col)
}


plotPhyloCoor_tmp <-
    function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL,
              direction = "rightwards", tip.height = NULL, ...)
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1)
        stop("found only one tip in the tree!")
    Nedge <- dim(x$edge)[1]
    if (any(tabulate(x$edge[, 1]) == 1))
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
    Nnode <- x$Nnode
    if (is.null(x$edge.length)) use.edge.length <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    if (phyloORclado) {
        ## changed by KS:
        yy <- numeric(Ntip + Nnode)
        x <- reorder(x)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        if (!is.null(tip.height)) {
            yy[TIPS] <- tip.height
        } 
        else yy[TIPS] <- 1:Ntip
    }
    
    xe <- x$edge
    ## first reorder the tree in cladewise order to avoid cophyloplot() hanging:
    ## x <- reorder(reorder(x), order = "pruningwise") ... maybe not needed anymore (EP)
    x <- reorder(x, order = "postorder")
    ereorder <- match(x$edge[, 2], xe[, 2])
    
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length)
                node.pos <- 2
    }
        if (node.pos == 1)
            yy <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
                     as.integer(x$edge[, 1]), as.integer(x$edge[,
                                                            2]), as.integer(Nedge), as.double(yy),
                     PACKAGE = "ape")[[6]]
        else {
            ans <- .C("node_height_clado", as.integer(Ntip),
                      as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
                                                                                    2]), as.integer(Nedge), double(Ntip + Nnode),
                      as.double(yy), PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2)
                xx <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                         as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                         as.integer(Nedge), double(Ntip + Nnode), 1L,
                         PACKAGE = "ape")[[6]] - 1
            xx <- max(xx) - xx
        } else {
            xx <- .C("node_depth_edgelength", as.integer(Ntip),
                     as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
                                                                                   2]), as.integer(Nedge), as.double(x$edge.length),
                     double(Ntip + Nnode), PACKAGE = "ape")[[7]]
        }
    }
    ##if (type == "fan") {
    ##    TIPS <- xe[which(xe[, 2] <= Ntip), 2]
    ##    xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
    ##    theta <- double(Ntip)
    ##    theta[TIPS] <- xx
    ##    theta <- c(theta, numeric(Nnode))
    ##    theta <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), theta, DUP = FALSE, PACKAGE = "ape")[[6]]
    ##    if (use.edge.length) {
    ##        r <- .C("node_depth_edgelength", as.integer(Ntip),
    ##            as.integer(Nnode), as.integer(x$edge[, 1]), as.integer(x$edge[,
    ##              2]), as.integer(Nedge), as.double(x$edge.length),
    ##            double(Ntip + Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    ##    }
    ##    else {
    ##        r <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
    ##            as.integer(x$edge[, 1]), as.integer(x$edge[,
    ##              2]), as.integer(Nedge), double(Ntip + Nnode),
    ##            DUP = FALSE, PACKAGE = "ape")[[6]]
    ##        r <- 1/r
    ##    }
    ##    xx <- r * cos(theta)
    ##    yy <- r * sin(theta)
    ##}
    ##if (type == "unrooted") {
    ##    XY <- if (use.edge.length)
    ##        unrooted.xy(Ntip, Nnode, x$edge, x$edge.length)
    ##    else unrooted.xy(Ntip, Nnode, x$edge, rep(1, Nedge))
    ##    xx <- XY$M[, 1] - min(XY$M[, 1])
    ##    yy <- XY$M[, 2] - min(XY$M[, 2])
    ##}
    ##if (type == "radial") {
    ##    X <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), double(Ntip + Nnode), DUP = FALSE,
    ##        PACKAGE = "ape")[[6]]
    ##    X[X == 1] <- 0
    ##    X <- 1 - X/Ntip
    ##    yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
    ##    Y <- .C("node_height", as.integer(Ntip), as.integer(Nnode),
    ##        as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
    ##        as.integer(Nedge), as.double(yy), DUP = FALSE, PACKAGE = "ape")[[6]]
    ##    xx <- X * cos(Y)
    ##    yy <- X * sin(Y)
    ##}
    if (phyloORclado && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    cbind(xx, yy)
}



#' Plots a densiTree.
#' 
#' An R function to plot trees similar to those produced by DensiTree.
#' 
#' If no consensus tree is provided \code{densiTree} computes a consensus tree, 
#' and if the input trees have different labels a mrp.supertree as a backbone. 
#' This should avoid too many unnecessary crossings of edges.  
#' Trees should be rooted, other wise the output may not be visually pleasing.
#' 
#' @param x an object of class \code{multiPhylo}.
#' @param type a character string specifying the type of phylogeny, so far
#' "cladogram" (default) or "phylogram" (the default) are supported.
#' @param alpha parameter for semi-transparent colors.
#' @param consensus A tree or character vector which is used to define the order 
#' of the tip labels.
#' @param direction a character string specifying the direction of the tree. 
#' Four values are possible: "rightwards" (the default), "leftwards", "upwards", 
#' and "downwards".
#' @param optim not yet used.
#' @param scaleX scale trees to have identical heights.
#' @param col a skalar or vector giving the colours used to draw the edges for 
#' each plotted phylogeny. These are taken to be in the same order than input 
#' trees x. If fewer colours are given than the number of trees, then the 
#' colours are recycled.
#' @param width edge width.
#' @param lty line type.
#' @param cex a numeric value giving the factor scaling of the tip labels.
#' @param font an integer specifying the type of font for the labels: 1 (plain text),
#'  2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param tip.color color of the tip labels. 
#' @param adj a numeric specifying the justification of the text strings of the 
#' labels: 0 (left-justification), 0.5 (centering), or 1 (right-justification). 
#' @param srt a numeric giving how much the labels are rotated in degrees.
#' @param underscore a logical specifying whether the underscores in tip labels 
#' should be written as spaces (the default) or left as are (if TRUE).  
#' @param label.offset a numeric giving the space between the nodes and the tips of the
#' phylogeny and their corresponding labels.
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
#'    function(x) upgma(dist.hamming(x)), bs=25)
#' # cladogram nice to show topological differences
#' densiTree(bs, type="cladogram", col="blue")
#' densiTree(bs, type="phylogram", col="green", direction="downwards", width=2)
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
densiTree <- function(x, type="cladogram", alpha=1/length(x), consensus=NULL, 
    direction="rightwards", optim=FALSE, scaleX=FALSE, col=1, width=1, lty=1,
    cex=.8, font=3, tip.color=1, adj=0, srt=0, underscore = FALSE, 
    label.offset=0, ...) {
  if(!inherits(x,"multiPhylo"))stop("x must be of class multiPhylo")

  if(is.character(consensus)){ 
      consensus <- stree(length(consensus), tip.label = consensus)
      consensus$edge.length <- rep(1.0, nrow(consensus$edge))
  }      
  if(is.null(consensus)){
      # unroot(midpoint(superTree(x)))
      consensus <- tryCatch(consensus(x, p=.5), error = function(e)unroot(midpoint(superTree(x))))   
  }      
  if(inherits(consensus,"multiPhylo")) consensus = consensus[[1]]


  type <- match.arg(type, c("phylogram", "cladogram"))
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards", 
    "downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")
  
  consensus = reorder(consensus)
  nTip = as.integer(length(consensus$tip.label))
#  ie <- match(1:nTip, consensus$edge[, 2])
  e2 <- consensus$edge[,2]
#  tiporder = e2[e2<=nTip]  
  consensus$tip.label <- consensus$tip.label[e2[e2<=nTip]]
  consensus$edge[e2<=nTip,2] <- as.integer(1L:nTip)
  consensus = reorder(consensus, "postorder")
  
#  ctmp <- consensus$edge[,2]
#  clabels <- consensus$tip.label[ ctmp[ctmp<cNtip] ] 

        
  x <- tryCatch(.compressTipLabel(x, ref = consensus$tip.label), error = function(e)x) #
  compressed <- ifelse(is.null(attr(x, "TipLabel")), FALSE, TRUE)
  
  maxBT = max(getAges(x))
  if(scaleX) maxBT=1.0
  label = rev(pretty(c(maxBT,0)))
  maxBT = max(label)
  xy = plotPhyloCoor_tmp(consensus, direction=direction, ...)
  yy = xy[,2]
  
  plot.new() 
  tl = which.max(nchar(consensus$tip.label))
  sw <- strwidth(consensus$tip.label[tl],cex=cex) * 1.1
   
  if(direction=="rightwards"){
    plot.window(xlim=c(0, 1.0+sw), ylim=c(0, nTip+1))
    axis(side=1,at=seq(0,1.0, length.out=length(label)), labels=label)
    #text(x=rep(1.0,Ntip(consensus)),y=xy[1:nTip, 2],labels=consensus$tip.label,pos=4,cex=cex, font=font)  
  }
  if(direction=="leftwards"){
      plot.window(xlim=c(0-sw, 1.0), ylim=c(0, nTip+1))
      axis(side=1,at=seq(0,1.0, length.out=length(label)), labels=rev(label))
      #text(x=rep(0,Ntip(consensus)),y=xy[1:nTip, 2],labels=consensus$tip.label,pos=2,cex=cex, font=font)  
  }
  if(direction=="downwards"){
      plot.window(xlim=c(0, nTip+1), ylim=c(0-sw, 1.0))
      axis(side=2,at=seq(0,1.0, length.out=length(label)), labels=rev(label))
      #text(x=xy[1:nTip,1],y=rep(0,Ntip(consensus)),labels=consensus$tip.label,pos=1,cex=cex, font=font, srt=90, adj=1)  
  }
  if(direction=="upwards"){
      plot.window(xlim=c(0, nTip+1), ylim=c(0, 1.0+sw))
      axis(side=2,at=seq(0,1.0, length.out=length(label)), labels=label)
      #text(x=xy[1:nTip,1],y=rep(1.0,Ntip(consensus)),labels=consensus$tip.label,pos=3,cex=cex, font=font, srt=90)  
  }
  if (is.expression(consensus$tip.label)) 
      underscore <- TRUE
  if (!underscore) 
      consensus$tip.label <- gsub("_", " ", consensus$tip.label)
 
  add_tiplabels(xy, consensus$tip.label, direction, adj=adj, font=font, srt=srt, 
                cex=cex, col=tip.color, label.offset=label.offset)
  
  col <- rep(col, length.out = length(x))
  tiporder <- NULL 
#  if(compressed) tiporder <- match(attr(x, "TipLabel"), consensus$tip.label)
#  tip.order = yy[1:nTip]
  for (treeindex in 1:length(x)) {
    tmp <- reorder(x[[treeindex]], "postorder")
    if(!compressed) tiporder <- match(tmp$tip.label, consensus$tip.label)
    xy <- plotPhyloCoor_tmp(tmp, tip.height=tiporder, direction=direction, ...)
    xx = xy[,1]
    yy = xy[,2]
    
    if(horizontal){  
      if(scaleX) xx <- xx/max(xx)
      else xx <- xx/maxBT 
      if(direction=="rightwards") xx <- xx + (1.0 - max(xx))
    }
    else{
      if(scaleX) yy <- yy/max(yy)
      else yy <- yy/maxBT 
      if(direction=="upwards")yy <- yy + (1.0 - max(yy))
    }
    e1 <- tmp$edge[,1]
    if(type=="cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color=
      adjustcolor(col[treeindex], alpha.f=alpha), edge.width=width, edge.lty=lty)
    if(type=="phylogram"){
      Ntip <- min(e1)-1L 
      Nnode <- tmp$Nnode 
      phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, horizontal, edge.color=
        adjustcolor(col[treeindex], alpha.f=alpha), edge.width=width, edge.lty=lty) 
    }
  }  
}



