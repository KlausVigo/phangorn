# some trigonemetric functions
rad2deg <- function(rad) (rad * 180) / (pi)
deg2rad <- function(deg) (deg * pi) / (180)

# circular mean
# https://en.wikipedia.org/wiki/Mean_of_circular_quantities
circ.mean <- function(deg) {
  rad.m <- (deg * pi) / (180)
  mean.cos <- mean(cos(rad.m))
  mean.sin <- mean(sin(rad.m))

  theta <- rad2deg(atan(mean.sin / mean.cos))
  if (mean.cos < 0) theta <- theta + 180
  if ((mean.sin < 0) & (mean.cos > 0)) theta <- theta + 360
  theta
}


spl2angle <- function(x) {
  l <- length(attr(x, "labels"))
  ord <- 1:l
  if (!is.null(attr(x, "cycle"))) ord <- attr(x, "cycle")
  x <- changeOrder(x, attr(x, "labels")[ord])
  y <- lapply(x, function(x, l) (x - 1) / l * 360, l = l)
  angle <- vapply(y, circ.mean, 0) |> deg2rad()
  # angle <- ((vapply(x, sum, 0) / lengths(x) - 1) / l ) * 2*pi
  # kreis2kart(attr(x, "weight"), angle)
  angle
}


coords.equal.angle <- function(obj) {
  if (is.null(attr(obj, "order")) || (attr(obj, "order") == "postorder"))
    obj <- reorder.networx(obj)
  spl <- obj$splits
  spl <- SHORTwise(spl) #, length(obj$tip.label))
  l <- length(obj$edge.length)
  #  ind1 <- which(!duplicated(obj$splitIndex))
  n <- max(obj$edge)
  angle <- spl2angle(spl)
  weight <- attr(spl, "weight")
  k <- matrix(0, max(obj$splitIndex), 2)

  res <- matrix(0, max(obj$edge), 2)
  for (i in 1:l) { # unique(obj$splitIndex)
    j <- obj$edge[i, 1]
    m <- obj$edge[i, 2]
    p <- obj$splitIndex[i]
    res[m, ] <- res[j, ] + kreis2kart(weight[p], angle[p])
  }
  res
}


#' @rdname phangorn-internal
#' @export
coords <- function(obj, dim = "3D") {
  #    if(is.null(attr(obj,"order")) || (attr(obj, "order")=="postorder") )
  #        obj = reorder.networx(obj)
  if (dim == "equal_angle") return(coords.equal.angle(obj))

  l <- length(obj$edge.length)
  ind1 <- which(!duplicated(obj$splitIndex))

  n <- max(obj$edge)
  adj <- spMatrix(n, n, i = obj$edge[, 2], j = obj$edge[, 1],
                  x = rep(1, length(obj$edge.length)))
  g <- graph_from_adjacency_matrix(adj, "undirected")
  if (dim == "3D") {
    coord <- layout_nicely(g, dim = 3)
    k <- matrix(0, max(obj$splitIndex), 3)
    for (i in ind1) {
      tmp <- coord[obj$edge[i, 2], ] - coord[obj$edge[i, 1], ]
      k[obj$splitIndex[i], ] <- kart2kugel(tmp[1], tmp[2], tmp[3])
    }
    k[obj$splitIndex[ind1], 1] <- obj$edge.length[ind1]

    res <- matrix(0, vcount(g), 3)
    for (i in 1:l) {
      j <- obj$edge[i, 1]
      m <- obj$edge[i, 2]
      p <- obj$splitIndex[i]
      res[m, ] <- res[j, ] + kugel2kart(k[p, 1], k[p, 2], k[p, 3])
    }
  }
  else {
    coord <- layout_nicely(g, dim = 2)
    k <- matrix(0, max(obj$splitIndex), 2)
    tmp <- coord[obj$edge[ind1, 2], ] - coord[obj$edge[ind1, 1], ]
    k[obj$splitIndex[ind1], ] <- kart2kreis(tmp[,1], tmp[,2])
    #    for (i in ind1) {
    #      tmp <- coord[obj$edge[i, 2], ] - coord[obj$edge[i, 1], ]
    #      k[obj$splitIndex[i], ] <- kart2kreis(tmp[1], tmp[2])
    #    }
    k[obj$splitIndex[ind1], 1] <- obj$edge.length[ind1]
    res <- matrix(0, vcount(g), 2)
    for (i in 1:l) {
      j <- obj$edge[i, 1]
      m <- obj$edge[i, 2]
      p <- obj$splitIndex[i]
      res[m, ] <- res[j, ] + kreis2kart(k[p, 1], k[p, 2])
    }
  }
  res
}


kart2kugel <- function(x, y, z) {
  r <- sqrt(x * x + y * y + z * z)
  alpha <- atan(sqrt(x * x + y * y) / z)
  if (z < 0) alpha <- alpha + pi
  beta <- atan(y / x)
  if (x < 0) beta <- beta + pi
  c(r, alpha, beta)
}


kart2kreis <- function(x, y) {
  r <- sqrt(x * x + y * y)
  alpha <- atan(y / x)
  #if (x < 0) alpha <- alpha + pi
  if (any(x < 0)) alpha[x < 0] <- alpha[x < 0] + pi
  cbind(r, alpha)
}


kreis2kart <- function(r, alpha) {
  c(r * cos(alpha), r * sin(alpha))
  #    if(length(r)>1) return(matrix(c(r*cos(alpha), r*sin(alpha)), ncol=2))
  # 	else return(c(r*cos(alpha), r*sin(alpha)))
}


kugel2kart <- function(r, alpha, beta) {
  x <- r * sin(alpha) * cos(beta)
  y <- r * sin(alpha) * sin(beta)
  z <- r * cos(alpha)
  c(x, y, z)
}


edgeLabels <- function(xx, yy, zz = NULL, edge) {
  XX <- (xx[edge[, 1]] + xx[edge[, 2]]) / 2
  YY <- (yy[edge[, 1]] + yy[edge[, 2]]) / 2
  if (!is.null(zz)) {
    ZZ <- (zz[edge[, 1]] + zz[edge[, 2]]) / 2
    return(cbind(XX, YY, ZZ))
  }
  cbind(XX, YY)
}


rotate_matrix <- function(x, theta){
  rot_matrix <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
                       2, 2, byrow = TRUE)
  x %*% rot_matrix
}


#' plot phylogenetic networks
#'
#' So far not all parameters behave the same on the the \code{rgl} \code{"3D"}
#' and basic graphic \code{"2D"} device.
#'
#' Often it is easier and safer to supply vectors of graphical parameters for
#' splits (e.g. splits.color) than for edges. These overwrite values edge.color.
#'
#' @param x an object of class \code{"networx"}
#' @param type "3D" to plot using rgl or "equal angle" and "2D" in the normal
#' device.
#' @param use.edge.length a logical indicating whether to use the edge weights
#' of the network to draw the branches (the default) or not.
#' @param show.tip.label a logical indicating whether to show the tip labels on
#' the graph (defaults to \code{TRUE}, i.e. the labels are shown).
#' @param show.edge.label a logical indicating whether to show the tip labels
#' on the graph.
#' @param edge.label an additional vector of edge labels (normally not needed).
#' @param show.node.label a logical indicating whether to show the node labels
#' (see example).
#' @param node.label an additional vector of node labels (normally not needed).
#' @param show.nodes a logical indicating whether to show the nodes (see
#' example).
#' @param tip.color the colors used for the tip labels.
#' @param edge.color the colors used to draw edges.
#' @param edge.width the width used to draw edges.
#' @param edge.lty a vector of line types.
#' @param split.color the colors used to draw edges.
#' @param split.width the width used to draw edges.
#' @param split.lty a vector of line types.
#' @param font an integer specifying the type of font for the labels: 1 (plain
#' text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#' @param cex a numeric value giving the factor scaling of the labels.
#' @param cex.node.label a numeric value giving the factor scaling of the node
#' labels.
#' @param cex.edge.label a numeric value giving the factor scaling of the edge
#' labels.
#' @param col.node.label the colors used for the node labels.
#' @param col.edge.label the colors used for the edge labels.
#' @param font.node.label the font used for the node labels.
#' @param font.edge.label the font used for the edge labels.
#' @param underscore a logical specifying whether the underscores in tip labels
#' should be written as spaces (the default) or left as are (if TRUE).
#' @param angle rotate the plot.
#' @param digits if edge labels are numerical a positive integer indicating how
#' many significant digits are to be used.
#' @param \dots Further arguments passed to or from other methods.
#' @returns \code{plot.networx} returns invisibly a list with paramters of the
#' plot.
#' @rdname plot.networx
#' @note The internal representation is likely to change.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{consensusNet}}, \code{\link{neighborNet}},
#' \code{\link{splitsNetwork}}, \code{\link{hadamard}},
#' \code{\link{distanceHadamard}}, \code{\link{as.networx}},
#' \code{\link[ape]{evonet}}, \code{\link[ape]{as.phylo}},
#' \code{\link{densiTree}}, \code{\link[ape]{nodelabels}}
#' @references Dress, A.W.M. and Huson, D.H. (2004) Constructing Splits Graphs
#' \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics
#' (TCBB)}, \bold{1(3)}, 109--115
#'
#' Schliep, K., Potts, A. J., Morrison, D. A. and Grimm, G. W. (2017),
#' Intertwining phylogenetic trees and networks. \emph{Methods Ecol Evol}.
#' \bold{8}, 1212--1220. doi:10.1111/2041-210X.12760
#' @keywords plot
#' @importFrom igraph make_graph
#' @examples
#'
#' set.seed(1)
#' tree1 <- rtree(20, rooted=FALSE)
#' sp <- as.splits(rNNI(tree1, n=10))
#' net <- as.networx(sp)
#' plot(net)
#' plot(net, direction="axial")
#' \dontrun{
#' # also see example in consensusNet
#' example(consensusNet)
#' }
#' @importFrom igraph graph_from_adjacency_matrix vcount topo_sort layout_nicely
#' @method plot networx
#' @export
plot.networx <- function(x, type = "equal angle", use.edge.length = TRUE,
                         show.tip.label = TRUE, show.edge.label = FALSE,
                         edge.label = NULL, show.node.label = FALSE,
                         node.label = NULL, show.nodes = FALSE,
                         tip.color = "black", edge.color = "black",
                         edge.width = 3, edge.lty = 1, split.color = NULL,
                         split.width = NULL, split.lty = NULL, font = 3,
                         cex = par("cex"), cex.node.label = cex,
                         cex.edge.label = cex, col.node.label = tip.color,
                         col.edge.label = tip.color, font.node.label = font,
                         font.edge.label = font, underscore = FALSE,
                         angle=0, digits=3, ...) {
  type <- match.arg(type, c("equal angle", "3D", "2D"))
  if (use.edge.length == FALSE){
    x$edge.length[] <- 1
    attr(x$splits, "weight") <- rep(1, length(x$splits))
  }
  nTips <- length(x$tip.label)
  conf <- attr(x$splits, "confidences")
  index <- x$splitIndex
  if(!is.null(edge.label) && is.numeric(edge.label)) edge.label <- prettyNum(edge.label)
  if (is.null(edge.label) && !is.null(conf)) {
    conf <- conf[index]
    if(is.numeric(conf)) conf <- prettyNum(format(conf, digits=digits))
    if (!is.null(x$translate)) conf[match(x$translate$node, x$edge[, 2])] <- ""
    else conf[x$edge[, 2] <= nTips] <- ""
    edge.label <- conf
  }
  if (is.null(node.label)) node.label <- as.character(1:max(x$edge))
  if (show.tip.label) node.label[1:nTips] <- ""
  if (show.tip.label){
    if (is.expression(x$tip.label)) underscore <- TRUE
    if (!underscore) x$tip.label <- gsub("_", " ", x$tip.label)
  }

  lspl <- max(x$splitIndex)
  if (!is.null(split.color)) {
    if (length(split.color) != lspl)
      stop("split.color must be same length as splits")
    else edge.color <- split.color[x$splitIndex]
  }
  if (!is.null(split.width)) {
    if (length(split.width) != lspl)
      stop("split.color must be same length as splits")
    else edge.width <- split.width[x$splitIndex]
  }
  if (!is.null(split.lty)) {
    if (length(split.lty) != lspl)
      stop("split.color must be same length as splits")
    else edge.lty <- split.lty[x$splitIndex]
  }

  chk <- FALSE

  if (type == "3D") chk <- requireNamespace("rgl", quietly = TRUE)
  if (!chk && type == "3D") {
    warning("type='3D' requires the package 'rgl', plotting in '2D' instead!\n")
    type <- "2D"
  }
  # use precomputed vertices when available
  coord <- NULL
  if (!is.null(x$.plot)) coord <- x$.plot$vertices

  if (type == "3D") {
    if (is.null(coord) || ncol(coord) != 3)
      coord <- coords(x, dim = "3D")
    plotRGL(coord, x, show.tip.label = show.tip.label,
            show.edge.label = show.edge.label, edge.label = edge.label,
            show.node.label = show.node.label, node.label = node.label,
            show.nodes = show.nodes, tip.color = tip.color,
            edge.color = edge.color, edge.width = edge.width, font = font,
            cex = cex, cex.node.label = cex.node.label,
            cex.edge.label = cex.edge.label, col.node.label = col.node.label,
            col.edge.label = col.edge.label, font.node.label = font.node.label,
            font.edge.label = font.edge.label)
  }
  else {
    if (is.null(coord) || ncol(coord) != 2) {
      if (type == "equal angle") coord <- coords.equal.angle(x)
      else coord <- coords(x, dim = "2D")
    }
    if(angle != 0){
      angle <- angle * pi/180 #
      coord <- rotate_matrix(coord, angle)
    }
    plot2D(coord, x, show.tip.label = show.tip.label,
           show.edge.label = show.edge.label, edge.label = edge.label,
           show.node.label = show.node.label, node.label = node.label,
           #show.nodes = show.nodes,
           tip.color = tip.color, edge.color = edge.color,
           edge.width = edge.width, edge.lty = edge.lty, font = font, cex = cex,
           cex.node.label = cex.node.label, cex.edge.label = cex.edge.label,
           col.node.label = col.node.label, col.edge.label = col.edge.label,
           font.node.label = font.node.label, font.edge.label = font.edge.label,
           add = FALSE, ...)
  }
  x$.plot <- list(vertices = coord, edge.color = edge.color,
                  edge.width = edge.width, edge.lty = edge.lty)
  invisible(x)
}


plotRGL <- function(coords, net, show.tip.label = TRUE, show.edge.label = FALSE,
                    edge.label = NULL, show.node.label = FALSE,
                    node.label = NULL, show.nodes = FALSE, tip.color = "blue",
                    edge.color = "grey", edge.width = 3, font = 3,
                    cex = par("cex"), cex.node.label = cex,
                    cex.edge.label = cex, col.node.label = tip.color,
                    col.edge.label = tip.color, font.node.label = font,
                    font.edge.label = font, ...) {
  open3d <- rgl::open3d
  segments3d  <- rgl::segments3d
  spheres3d <- rgl::spheres3d
  texts3d <- rgl::texts3d

  edge <- net$edge

  x <- coords[, 1]
  y <- coords[, 2]
  z <- coords[, 3]

  nTips <- length(net$tip.label)

  segments3d(x[t(edge)], y[t(edge)], z[t(edge)],
             col = rep(edge.color, each = 2), lwd = edge.width)
  radius <- 0
  if (show.nodes) {
    radius <- sqrt( (max(x) - min(x))^2 + (max(y) - min(y))^2 +
                      (max(z) - min(z))^2) / 200
    spheres3d(x[1:nTips], y[1:nTips], z[1:nTips], radius = 2 * radius,
              color = "cyan")
    spheres3d(x[-c(1:nTips)], y[-c(1:nTips)], z[-c(1:nTips)], radius = radius,
              color = "magenta")
  }
  if (show.tip.label) {
    if (is.null(net$translate))
      texts3d(x[1:nTips] + 2.05 * radius, y[1:nTips], z[1:nTips],
              net$tip.label, color = tip.color, cex = cex, font = font)
    else
      texts3d(x[net$translate$node] + 2.05 * radius, y[net$translate$node],
              z[net$translate$node], net$tip.label, color = tip.color,
              cex = cex, font = font)
  }
  if (show.edge.label) {
    ec <- edgeLabels(x, y, z, edge)
    if (is.null(edge.label)) edge.label <- net$splitIndex
    # else edge.label = net$splitIndex
    texts3d(ec[, 1], ec[, 2], ec[, 3], edge.label, color = col.edge.label,
            cex = cex.edge.label, font = font.edge.label)
  }
  if (show.node.label) {
    texts3d(x, y, z, node.label, color = col.node.label, cex = cex.node.label,
            font = font.node.label)
  }
}


plot2D <- function(coords, net, show.tip.label = TRUE, show.edge.label = FALSE,
                   edge.label = NULL, show.node.label = FALSE,
                   node.label = NULL, tip.color = "blue", edge.color = "grey",
                   edge.width = 3, edge.lty = 1, font = 3, cex = par("cex"),
                   cex.node.label = cex,  cex.edge.label = cex,
                   col.node.label = tip.color, col.edge.label = tip.color,
                   font.node.label = font, font.edge.label = font,
                   add = FALSE, direction="horizontal", xlim=NULL, ylim=NULL,
                   ...) {
  direction <- match.arg(direction, c("horizontal", "axial"))
  edge <- net$edge
  label <- net$tip.label
  xx <- coords[, 1]
  yy <- coords[, 2]
  nTips <- length(label)
  offset <- max(nchar(label)) * 0.018 * cex * diff(range(xx))
  if(is.null(xlim)){
    xlim <- range(xx)
#    offset <- max(nchar(label)) * 0.018 * cex * diff(xlim)
    if (show.tip.label) xlim <- c(xlim[1] - offset, xlim[2] + offset)
  }
  if(is.null(ylim)){
    ylim <- range(yy)
    if (show.tip.label){
      if(direction=="axial"){
#        offset <- max(nchar(label)) * 0.018 * cex * diff(ylim)
        ylim <- c(ylim[1] - offset, ylim[2] + offset)
      } else ylim <- c(ylim[1] - 0.03 * cex * diff(ylim),
                     ylim[2] + 0.03 * cex * diff(ylim))
    }
  }
  if (!add) {
#    plot.new()
#    plot.window(xlim, ylim, asp = 1, ...)
    plot.default(0, type = "n", xlim = xlim, ylim = ylim, xlab = "",
                 ylab = "", axes = FALSE, asp = 1, ...)
  }
  cladogram.plot(edge, xx, yy, edge.color, edge.width, edge.lty)
  if (show.tip.label) {
    if (is.null(net$translate)) ind <- match(1:nTips, edge[, 2])
    else ind <- match(net$translate$node, edge[, 2])
    if(direction=="horizontal"){
      pos <- rep(4, nTips)
      XX <- xx[edge[ind, 1]] - xx[edge[ind, 2]]
      pos[XX > 0] <- 2
      YY <- yy[edge[ind, 1]] - yy[edge[ind, 2]]
      pos2 <- rep(3, nTips)
      pos2[YY > 0] <- 1
      # needed if tiplabels are not at internal nodes
      XX[is.na(XX)] <- 0
      YY[is.na(YY)] <- 0
      pos[abs(YY) > abs(XX)] <- pos2[abs(YY) > abs(XX)]
      if (is.null(net$translate)) text(xx[1:nTips], yy[1:nTips], labels = label,
                                       pos = pos, col = tip.color, cex = cex, font = font)
      else text(xx[net$translate$node], yy[net$translate$node], labels = label,
                pos = pos, col = tip.color, cex = cex, font = font)
    }
    else {
      XX <- xx[edge[ind, 2]] - xx[edge[ind, 1]]
      YY <- yy[edge[ind, 2]] - yy[edge[ind, 1]]
      angle <- kart2kreis(XX, YY)[,2]
      adj <- abs(angle) > pi/2
      angle <- angle * 180/pi # switch to degrees
      angle[adj] <- angle[adj] - 180
      adj <- as.numeric(adj)
      ## `srt' takes only a single value, so can't vectorize this:
      ## (and need to 'elongate' these vectors:)
      font <- rep(font, length.out = nTips)
      tip.color <- rep(tip.color, length.out = nTips)
      cex <- rep(cex, length.out = nTips)
      for (i in seq_along(label))
        text(xx[i], yy[i], label[i], font = font[i],
             cex = cex[i], srt = angle[i], adj = adj[i],
             col = tip.color[i])
    }
  }

  if (show.edge.label) {
    ec <- edgeLabels(xx, yy, edge = edge)
    if (is.null(edge.label)) edge.label <- net$splitIndex

    # show only one edge label
    em <- apply(ec, 1, function(x) max(abs(x)))
    si <- net$splitIndex
    for (i in unique(si)) {
      tmp <- si == i
      if (sum(tmp) > 1) {
        w <- which(tmp)
        wm <- which.max(em[w])
        edge.label[w[-wm]] <- ""
      }
    }

    text(ec[, 1], ec[, 2], labels = edge.label, col = col.edge.label,
         cex = cex.edge.label, font = font.edge.label)
  }
  if (show.node.label) {
    text(xx, yy, labels = node.label, col = col.node.label,
         cex = cex.node.label, font = font.node.label)
  }
  PP <- list(Ntip = nTips, type = "networx", edge = net$edge, xx = coords[, 1],
             yy = coords[, 2], x.lim=xlim, y.lim=ylim, align.tip.label=FALSE)
  assign("last_plot.phylo", PP, envir = .PlotPhyloEnv)
}


closest.edge <- function(x, y, P1, P2) {
  x1 <- P1[, 1]
  x2 <- P2[, 1]
  y1 <- P1[, 2]
  y2 <- P2[, 2]

  A <- sqrt( (x2 - x)^2 + (y2 - y)^2)    # d_BC
  B <- sqrt( (x1 - x)^2 + (y1 - y)^2)    # d_AC
  C <- sqrt( (x1 - x2)^2 + (y1 - y2)^2)  # d_AB
  # Kosinussatz
  alpha <- acos( (B^2 + C^2 - A^2) / (2 * B * C))
  beta <- acos( (A^2 + C^2 - B^2) / (2 * A * C))

  d <- abs( (y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    sqrt( (y2 - y1)^2 + (x2 - x1)^2)
  d[alpha > (pi / 2)] <- B[alpha > (pi / 2)]
  d[beta > (pi / 2)] <- A[beta > (pi / 2)]
  d
}

closest.node <- function(x, y, P) {
  x1 <- P[, 1]
  y1 <- P[, 2]
  d <- sqrt((x1 - x)^2 + (y1 - y)^2)
  d
}


#' Identify splits in a network
#'
#' \code{identify.networx} reads the position of the graphics pointer when the
#' mouse button is pressed. It then returns the split belonging to the edge
#' closest to the pointer. The network must be plotted beforehand.
#'
#' @param x an object of class \code{networx}
#' @param quiet a logical controlling whether to print a message inviting the
#' user to click on the tree.
#' @param \dots further arguments to be passed to or from other methods.
#' @return \code{identify.networx} returns a splits object.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[phangorn]{plot.networx}},
#'   \code{\link[graphics]{identify}}
#' @examples
#' \dontrun{
#' data(yeast)
#' dm <- dist.ml(yeast)
#' nnet <- neighborNet(dm)
#' plot(nnet)
#' identify(nnet) # click close to an edge
#' }
#' @importFrom graphics identify
#' @method identify networx
#' @export
identify.networx <- function(x, quiet = FALSE, ...) {
  if (!quiet)
    cat("Click close to a node or edge of the tree...\n")
  xy <- locator(1)
  if (is.null(xy))
    return(NULL)
  if (is.null(x$.plot)) {
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    edge <- lastPP$edge
    xx <- lastPP$xx
    yy <- lastPP$yy
    vertices <- cbind(xx, yy)
  }
  else {
    lastPP <- x$.plot
    edge <- x$edge
    vertices <- lastPP$vertices
  }
  P1 <- vertices[edge[, 1], , drop = FALSE]
  P2 <- vertices[edge[, 2], , drop = FALSE]
  d <- closest.edge(xy$x, xy$y, P1, P2)
  split <- x$splitIndex[which.min(d)]
  x$splits[split]
}
