#' Bootstrap
#'
#' \code{bootstrap.pml} performs (non-parametric) bootstrap analysis and
#' \code{bootstrap.phyDat} produces a list of bootstrapped data sets.
#' \code{plotBS} plots a phylogenetic tree with the bootstrap values assigned
#' to the (internal) edges.
#'
#' It is possible that the bootstrap is performed in parallel, with help of the
#' multicore package. Unfortunately the multicore package does not work under
#' windows or with GUI interfaces ("aqua" on a mac). However it will speed up
#' nicely from the command line ("X11").
#'
#' @param x an object of class \code{pml} or \code{phyDat}.
#' @param bs number of bootstrap samples.
#' @param trees return trees only (default) or whole \code{pml} objects.
#' @param multicore logical, whether models should estimated in parallel.
#' @param mc.cores The number of cores to use during bootstrap. Only supported
#' on UNIX-alike systems.
#' @param jumble logical, jumble the order of the sequences.
#' @param tip.dates	 A named vector of sampling times associated to the
#' tips/sequences. Leave empty if not estimating tip dated phylogenies.
#' @param \dots further parameters used by \code{optim.pml} or
#' \code{plot.phylo}.
#' @param FUN the function to estimate the trees.
#' @return \code{bootstrap.pml} returns an object of class \code{multi.phylo}
#' or a list where each element is an object of class \code{pml}. \code{plotBS}
#' returns silently a tree, i.e. an object of class \code{phylo} with the
#' bootstrap values as node labels. The argument \code{BStrees} is optional and
#' if not supplied the tree with labels supplied in the \code{node.label} slot.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{optim.pml}}, \code{\link{pml}},
#' \code{\link{plot.phylo}}, \code{\link{maxCladeCred}}
#' \code{\link{nodelabels}},\code{\link{consensusNet}} and
#' \code{\link{SOWH.test}} for parametric bootstrap
#' @references Felsenstein J. (1985) Confidence limits on phylogenies. An
#' approach using the bootstrap. \emph{Evolution} \bold{39}, 783--791
#'
#' Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe, M. D.,
#' De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s phylogenetic
#' bootstrap in the era of big data. \emph{Nature}, \bold{556(7702)}, 452--456.
#'
#' Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree
#' construction. \emph{Cladistics} \bold{1}, 266--278
#'
#' Penny D. and Hendy M.D. (1986) Estimating the reliability of evolutionary
#' trees. \emph{Molecular Biology and Evolution} \bold{3}, 403--417
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(Laurasiatherian)
#' dm <- dist.hamming(Laurasiatherian)
#' tree <- NJ(dm)
#' # NJ
#' set.seed(123)
#' NJtrees <- bootstrap.phyDat(Laurasiatherian,
#'      FUN=function(x)NJ(dist.hamming(x)), bs=100)
#' treeNJ <- plotBS(tree, NJtrees, "phylogram")
#'
#' # Maximum likelihood
#' fit <- pml(tree, Laurasiatherian)
#' fit <- optim.pml(fit, rearrangement="NNI")
#' set.seed(123)
#' bs <- bootstrap.pml(fit, bs=100, optNni=TRUE)
#' treeBS <- plotBS(fit$tree,bs)
#'
#' # Maximum parsimony
#' treeMP <- pratchet(Laurasiatherian)
#' treeMP <- acctran(treeMP, Laurasiatherian)
#' set.seed(123)
#' BStrees <- bootstrap.phyDat(Laurasiatherian, pratchet, bs = 100)
#' treeMP <- plotBS(treeMP, BStrees, "phylogram")
#' add.scale.bar()
#'
#' # export tree with bootstrap values as node labels
#' # write.tree(treeBS)
#' }
#'
#' @rdname bootstrap.pml
#' @export
bootstrap.pml <- function(x, bs = 100, trees = TRUE, multicore = FALSE,
                          mc.cores = NULL, tip.dates=NULL, ...) {
  if(.Platform$OS.type=="windows") multicore <- FALSE
  if (multicore && is.null(mc.cores)) mc.cores <- min(detectCores()-1L, 4L)
  if(multicore && mc.cores < 2L) multicore <- FALSE
  if(is.rooted(x$tree)){
    if(is.ultrametric(x$tree)) method <- "ultrametric"
    else method <- "tipdated"
  }
  else method <- "unrooted"

  extras <- match.call(expand.dots = FALSE)$...
  rearr <- c("optNni", "rearrangement")
  tmp <- pmatch(names(extras), rearr)
  tmp <- tmp[!is.na(tmp)]
  do_rearr <- FALSE
  if(length(tmp)>0){
    if(tmp==1){
      do_rearr <- extras$optNni
      if(is.name(do_rearr)) do_rearr <- as.logical(as.character(do_rearr))
    }
    if(tmp==2) do_rearr <- extras$rearrangement %in% c("NNI", "stochastic",
                                                       "ratchet")
  }
  is_ultrametric <- FALSE
  tmp <- pmatch("optRooted", names(extras))
  if(!is.na(tmp)){
    is_ultrametric <- extras$optRooted
  }
  data <- x$data
  weight <- attr(data, "weight")
  v <- rep(seq_along(weight), weight)
  ntips <- Ntip(x$tree)
  BS <- vector("list", bs)
  for (i in 1:bs) BS[[i]] <- tabulate(
      sample(v, replace = TRUE),
      length(weight)
    )
  pmlPar <- function(weights, fit, trees = TRUE, do_rearr, ...) {
    data <- fit$data
    tree <- fit$tree
    ind <- which(weights > 0)
    data <- getRows(data, ind)
    attr(data, "weight") <- weights[ind]
    fit <- update(fit, data = data)
    if(do_rearr){
      tree <- candidate_tree(data, method=method, bf=fit$bf, Q=fit$Q,
                             tip.dates = tip.dates)
      fit <- update(fit, tree = tree)
    }
    fit <- optim.pml(fit, ...)
    if (trees) {
      tree <- fit$tree
      return(tree)
    }
    attr(fit, "data") <- NULL
    fit
  }
  eval.success <- FALSE
  if(method=="tipdated") do_rearr <- FALSE
  if (!eval.success & multicore) {
    res <- mclapply(BS, pmlPar, x, trees = trees, do_rearr = do_rearr, ...,
                    mc.cores = mc.cores)
    eval.success <- TRUE
  }
  if (!eval.success) res <- lapply(BS, pmlPar, x, trees = trees,
                                   do_rearr = do_rearr, ...)
  if (trees) {
    class(res) <- "multiPhylo"
    res <- .compressTipLabel(res) # save memory
  }
  res
}

#' @rdname bootstrap.pml
#' @export
bootstrap.phyDat <- function(x, FUN, bs = 100, multicore = FALSE,
                             mc.cores = NULL, jumble = TRUE, ...) {
  if(.Platform$OS.type=="windows") multicore <- FALSE
  if (multicore && is.null(mc.cores)) mc.cores <- detectCores()
  weight <- attr(x, "weight")
  v <- rep(seq_along(weight), weight)
  BS <- vector("list", bs)
  for (i in 1:bs) BS[[i]] <- tabulate(sample(v, replace = TRUE), length(weight))
  if (jumble) {
    J <- vector("list", bs)
    l <- length(x)
    for (i in 1:bs) J[[i]] <- list(BS[[i]], sample(l))
  }
  fitPar <- function(weights, data, ...) {
    ind <- which(weights > 0)
    data <- getRows(data, ind)
    attr(data, "weight") <- weights[ind]
    FUN(data, ...)
  }
  fitParJumble <- function(J, data, ...) {
    ind <- which(J[[1]] > 0)
    data <- getRows(data, ind)
    attr(data, "weight") <- J[[1]][ind]
    data <- subset(data, J[[2]])
    FUN(data, ...)
  }
  if (multicore) {
    if (jumble) {
      res <- mclapply(J, fitParJumble, x, ..., mc.cores = mc.cores)
    } else {
      res <- mclapply(BS, fitPar, x, ..., mc.cores = mc.cores)
    }
  }
  else {
    if (jumble) {
      res <- lapply(J, fitParJumble, x, ...)
    } else {
      res <- lapply(BS, fitPar, x, ...)
    }
  }
  if (inherits(res[[1]], "phylo")) {
    class(res) <- "multiPhylo"
    res <- .compressTipLabel(res) # save memory
  }
  res
}


checkLabels <- function(tree, tip) {
  ind <- match(tree$tip.label, tip)
  if (any(is.na(ind)) | length(tree$tip.label) != length(tip)) {
    stop("tree has different labels")
  }
  tree$tip.label <- tip #tree$tip.label[ind]
  ind2 <- tree$edge[, 2] <= Ntip(tree)
  tree$edge[ind2, 2] <- ind[tree$edge[ind2, 2]]
  tree
}


#' Plotting trees with bootstrap values
#'
#' \code{plotBS} plots a phylogenetic tree with the bootstrap values assigned
#' to the (internal) edges. It can also used to assign bootstrap values to a
#' phylogenetic tree.
#'
#' \code{plotBS} can either assign the classical Felsenstein’s bootstrap
#' proportions (FBP) (Felsenstein (1985), Hendy & Penny (1985))  or the
#' transfer bootstrap expectation (TBE) of Lemoine et al. (2018). Using the
#' option \code{type=="n"} just assigns the bootstrap values and return the tree
#' without plotting it.
#'
#' @param tree The tree on which edges the bootstrap values are plotted.
#' @param BStrees a list of trees (object of class "multiPhylo").
#' @param type the type of tree to plot, one of "phylogram", "cladogram", "fan",
#' "unrooted", "radial" or "none". If type is "none" the tree is returned with
#' the bootstrap values assigned to the node labels.
#' @param method either "FBP" the classical bootstrap (default) or "TBE"
#' (transfer bootstrap)
#' @param bs.col color of bootstrap support labels.
#' @param bs.adj one or two numeric values specifying the horizontal and
#' vertical justification of the bootstrap labels.
#' @param digits integer indicating the number of decimal places.
#' @param p only plot support values higher than this percentage number
#' (default is 0).
#' @param \dots further parameters used by \code{plot.phylo}.
#' @param frame a character string specifying the kind of frame to be printed
#' around the bootstrap values. This must be one of "none" (the default),
#' "rect" or "circle".
#' @return \code{plotBS} returns silently a tree, i.e. an object of class
#' \code{phylo} with the bootstrap values as node labels. The argument
#' \code{BStrees} is optional and if not supplied the labels supplied
#' in the \code{node.label} slot will be used.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso  \code{\link{transferBootstrap}}, \code{\link{plot.phylo}},
#' \code{\link{maxCladeCred}}, \code{\link{nodelabels}},
#' \code{\link{consensus}}, \code{\link{consensusNet}}
#' @references Felsenstein J. (1985) Confidence limits on phylogenies. An
#' approach using the bootstrap. \emph{Evolution} \bold{39}, 783--791
#'
#' Lemoine, F., Entfellner, J. B. D., Wilkinson, E., Correia, D., Felipe, M. D.,
#' De Oliveira, T., & Gascuel, O. (2018). Renewing Felsenstein’s phylogenetic
#' bootstrap in the era of big data. \emph{Nature}, \bold{556(7702)}, 452--456.
#'
#' Penny D. and Hendy M.D. (1985) Testing methods evolutionary tree
#' construction. \emph{Cladistics} \bold{1}, 266--278
#'
#' Penny D. and Hendy M.D. (1986) Estimating the reliability of evolutionary
#' trees. \emph{Molecular Biology and Evolution} \bold{3}, 403--417
#' @examples
#' fdir <- system.file("extdata/trees", package = "phangorn")
#' # RAxML best-known tree with bipartition support (from previous analysis)
#' raxml.tree <- read.tree(file.path(fdir,"RAxML_bipartitions.woodmouse"))
#' # RAxML bootstrap trees (from previous analysis)
#' raxml.bootstrap <- read.tree(file.path(fdir,"RAxML_bootstrap.woodmouse"))
#' par(mfrow=c(1,2))
#' plotBS(raxml.tree,  raxml.bootstrap, "p")
#' plotBS(raxml.tree,  raxml.bootstrap, "p", "TBE")
#' @export
plotBS <- function(tree, BStrees, type = "phylogram",
                   method="FBP", bs.col = "black",
                   bs.adj = NULL, digits=3, p = 0, frame = "none", ...) {
  type <- match.arg(type, c("phylogram", "cladogram", "fan", "unrooted",
                            "radial", "none"))
  method <- match.arg(method, c("FBP", "TBE"))
  if (hasArg(BStrees)) {
    if(method=="FBP"){
      BStrees <- .uncompressTipLabel(BStrees) # check if needed
      if (any(is.rooted(BStrees))) BStrees <- unroot(BStrees)
      x <- prop.clades(tree, BStrees)
      x <- (x / length(BStrees)) * 100
      tree$node.label <- x
    }
    else {
      tree <- transferBootstrap(tree, BStrees)
      x <- tree$node.label
    }
  }
  else {
    if (is.null(tree$node.label)) stop("You need to supply 'trees' or the tree needs support-values as node.label")
    x <- tree$node.label
  }
  if(type=="none") return( tree )
    plot(tree, type = type, ...)

  label <- c(rep(0, length(tree$tip.label)), x)
  ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[ ,2 ]
  if (type == "phylogram" | type == "cladogram") {
    root <- getRoot(tree)
    label <- c(rep(0, length(tree$tip.label)), x)
    label[root] <- 0
    ind <- which(label > p)
    if (is.null(bs.adj)) {
      bs.adj <- c(1, 1)
    }
    if (length(ind) > 0) {
      if(is.numeric(label)) label <- round(label, digits = digits)
      nodelabels(
        text = label[ind], node = ind,
        frame = frame, col = bs.col, adj = bs.adj, ...
      )
    }
  }
  else {
    if (is.null(bs.adj)) {
      bs.adj <- c(0.5, 0.5)
    }
    ind2 <- which(label[ind] > p)
    if (length(ind2 > 0)) {
      if(is.numeric(label)) label <- round(label, digits = digits)
      edgelabels(label[ind][ind2], ind2,
        frame = frame,
        col = bs.col, adj = bs.adj, ...
      )
    }
  }
  invisible(tree)
}


cladeMatrix <- function(x, rooted = FALSE) {
  if (!rooted) x <- unroot(x)
  pp <- prop.part(x)
  pplabel <- attr(pp, "labels")
  if (!rooted) pp <- SHORTwise(pp)
  x <- .uncompressTipLabel(x)
  nnodes <- Nnode(x)
  class(x) <- NULL
  # nnodes <- sapply(x, Nnode)
  l <- length(x)
  from <- cumsum(c(1, nnodes[-l]))
  to <- cumsum(nnodes)
  ivec <- integer(to[l])
  pvec <- c(0, to)
  res <- vector("list", l)
  k <- 1
  for (i in 1:l) {
    ppi <- prop.part(x[[i]])
    if (!rooted) ppi <- SHORTwise(ppi)
    indi <- sort(fmatch(ppi, pp))
    ivec[from[i]:to[i]] <- indi
  }
  X <- sparseMatrix(i = ivec, p = pvec, dims = c(length(pp), l))
  list(X = X, prop.part = pp)
}


moving_average <- function(obj, window = 50) {
  fun <- function(x) {
    cx <- c(0, cumsum(x))
    (cx[(window + 1):length(cx)] - cx[1:(length(cx) - window)]) / (window)
  }
  res <- apply(obj$X, 1, fun)
  rownames(res) <- c()
}


