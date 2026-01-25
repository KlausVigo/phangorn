#' @rdname pml
#' @export
logLik.pml <- function(object, ...) {
  res <- object$logLik
  attr(res, "df") <- object$df
  class(res) <- "logLik"
  res
}


#' @export
AICc <- function(object, ...)
  UseMethod("AICc")


#' @export
AICc.pml <- function(object, ...) {
  n <- sum(object$weight)
  k <- object$df
  if (k >= (n - 1)) return(NaN)
  res <- AIC(object)
  res +   (2 * k * (k + 1)) / (n - k - 1)
}


#' @export
BIC.pml <- function(object, ...) {
  res <- AIC(object, k = log(sum(object$weight)))
  res
}

#' @rdname pml
#' @export
anova.pml <- function(object, ...) {
  X <- c(list(object), list(...))
  fun <- function(x) {
    tmp <- logLik(x)
    c(tmp[1], attr(tmp, "df"))
  }
  DF <- t(sapply(X, fun))
  dev <- c(NA, 2 * diff(DF[, 1]))
  ddf <- c(NA, diff(DF[, 2]))
  table <- data.frame(DF, ddf, dev, pchisq(dev, ddf, lower.tail = FALSE))
  dimnames(table) <- list(seq_along(X), c("Log lik.", "Df", "Df change",
                                          "Diff log lik.", "Pr(>|Chi|)"))
  structure(table, heading = "Likelihood Ratio Test Table",
            class = c("anova", "data.frame"))
}


#' @rdname pml
#' @export
vcov.pml <- function(object, ...) {
  FI <- score(object, FALSE)[[2]]
  l <- dim(FI)[1]
  res <- try(solve(FI))
  if (inherits(res, "try-error")) {
    cat("Covariance is ill-conditioned !! \n")
    res <- solve(FI + diag(l) * 1e-8)
  }
  res
}

#' @rdname pml
#' @export
print.pml <- function(x, ...) {
  model <- guess_model(x)
  cat("model:", model, "\n")
  cat("loglikelihood:", x$logLik, "\n")
  w <- x$weight
  w <- w[w > 0]
  type <- attr(x$data, "type")
  levels <- attr(x$data, "levels")
  nc <- attr(x$data, "nc")
  ll0 <- sum(w * log(w / sum(w)))
  cat("unconstrained loglikelihood:", ll0, "\n")
  cat("Total tree length:",
      sum(x$tree$edge.length) * x$rate, "\n\t(expected number of substituions per site)\n")
  cat("Minimal tree length:",
      parsimony(x$tree, x$data) / sum(attr(x$data, "weight")), "\n\t(observed substitutions per site)\n")
  if (x$inv > 0) cat("Proportion of invariant sites:", x$inv, "\n")
  if (x$k > 1) {
    cat("Model of rate heterogeneity: ")
    if(x$site.rate=="gamma") cat("Discrete gamma model\n")
    if(x$site.rate=="free_rate") cat("Free rate model\n")
    if(x$site.rate=="gamma_quadrature") cat("Discrete gamma model (quadrature) \n")
    if(x$site.rate=="gamma_phangorn") cat("Discrete gamma model (phangorn) \n")
    cat("Number of rate categories:", x$k, "\n")
    if(x$site.rate!="free_rate") cat("Shape parameter:", x$shape, "\n")
    rate <- x$g
    prop <- x$w
    if (x$inv > 0) {
      rate <- c(0, rate)
      prop <- c(x$inv, prop)
    }
    rw <- cbind(Rate=rate, Proportion=prop)
    row.names(rw) <- as.character(seq_len(nrow(rw)))
    print(rw)
  }
  if(!is.null(x$method) && x$method == "tipdated") cat("\nRate:", x$rate, "\n")
  if (type == "AA") cat("Rate matrix:", x$model, "\n")
  if (type == "DNA") {

    lev <- attr(x$data, "levels")
    tri_ind <- function (n.rows)
    {
      seqi <- seq.int(n.rows - 1L)
      hi <- sequence(rev(seqi), from = seqi + 1)
      lo <- rep.int(seqi, rev(seqi))
      cbind(lo, hi, deparse.level = 0)
    }
    IND <- tri_ind(length(lev))
    cat("\nRates:\n")
    for(i in seq_len(nrow(IND))){
      cat(lev[IND[i,1]], "<->", lev[IND[i,2]], ":", x$Q[i], "\n")
    }
#    cat("\nRate matrix:\n")
#    QM <- matrix(0, nc, nc, dimnames = list(levels, levels))
#    QM[lower.tri(QM)] <- x$Q
#    QM <- QM + t(QM)
#    print(QM)
    cat("\nBase frequencies:  \n")
    bf <- x$bf
    names(bf) <- levels
    print(bf) #cat(bf, "\n")
  }
  if (type == "CODON") {
    cat("dn/ds:", x$dnds, "\n")
    cat("ts/tv:", x$tstv, "\n")
    cat("Freq:", x$frequencies, "\n")
  }
  if (type == "USER" & length(x$bf) < 11) {
    cat("\nRate matrix:\n")
    QM <- matrix(0, nc, nc, dimnames = list(levels, levels))
    QM[lower.tri(QM)] <- x$Q
    QM <- QM + t(QM)
    print(QM)
    cat("\nBase frequencies:  \n")
    bf <- x$bf
    names(bf) <- levels
    print(bf) #cat(bf, "\n")
  }
  if (!isTRUE(all.equal(x$rate, 1))) cat("\nRate:", x$rate, "\n")
  invisible(x)
}


#' Export pml objects
#'
#' \code{write.pml} writes out the ML tree and the model parameters.
#'
#' \code{write.pml} creates several files. It exports the alignment as fasta
#' file. It writes out the ML tree in a newick file and the estimates parameters
#' in a txt file. It should be possible to (re-)create the pml object up to
#' numerical inaccuracies and this is possible with the *.rds file.
#' If bootstrap trees exist these are additionally exported in a compressed
#' nexus file.
#' Additionally several plots are returned. The maximum likelihood tree, with
#' support values, if these are available. If an bootstrapped trees exist, a
#' consensus tree, a consensus network (< 200 tips) and terrace plot.
#' And last but not least the distribution of the rates. It might be better to
#' adopt these on the dataset.
#'
#' @param x an object of class ancestral.
#' @param file a file name. File endings are added.
#' @param save_rds logical, if TRUE saves the pml object as a rds file,
#' otherwise the alignment is saved as a fasta file.
#' @param digits default is 10, i.e. edge length for the bootstrap trees are
#' exported. For digits larger smaller than zero no edge length are exported.
## @param chi_sq logical, if TRUE performs $Chi^2$-test to check if sequences
## have similar state composition.
#' @param ... Further arguments passed to or from other methods.
#' @returns \code{write.pml}  returns the input x invisibly.
#' @seealso \code{\link{ancestral.pml}}, \code{\link{plotAnc}}
#' @examples
#' data(woodmouse)
#' fit <- pml_bb(woodmouse, "JC", rearrangement = "none")
#' write.pml(fit, "woodmouse")
#' unlink(c("woodmouse.txt", "woodmouse_tree.nwk", "woodmouse_align.fasta"))
#' @importFrom utils citation
#' @export
write.pml <- function(x, file="pml", save_rds=TRUE, digits=10, ...){
#  digits <- -1
#  if (hasArg("digits")) digits <- list(...)$digits
  write.tree(x$tree, file=paste0(file, "_tree.nwk"))
  if(save_rds) saveRDS(x, file=paste0(file, ".rds"))
  write.phyDat(x$data, file=paste0(file, "_align.fasta"), format="fasta")

  ntips <- Ntip(x$tree)
  my_heigtht <- max(7, ntips*50 / 300 )
  if(!is.null(x$bs)){
    zzfil <- paste0(file, "_bs.nex.gz")
    zz <- gzfile(zzfil, "w")
    write.nexus(x$bs, file=zz, digits=digits)
    close(zz)
    if(ntips < 3000){
      pdf(file=paste0(file, "_consensus.pdf"), height=my_heigtht)
      ctree <- consensus(x$bs, p=0.5)
      par(mar=c(1,1,1,1))
      plot(ctree, main="Majority consensus tree", node.depth = 2)
      dev.off()
    }
    if(ntips < 200){
      pdf(file=paste0(file, "_consensusNet.pdf"), height=my_heigtht)
      ctree <- consensus(x$bs, p=0.5)
      cnet <- consensusNet(x$bs, prob=0.3)
      plot(cnet, main="Consensus network", direction="axial")
      dev.off()
    }
    pdf(file=paste0(file, "_terraces.pdf"), height=my_heigtht)
    terraces(fit, pkg="scatterplot3d")
    dev.off()
  }
  if(ntips < 3000){
    pdf(file=paste0(file, "_tree.pdf"), height=my_heigtht)
    par(mar=c(1,1,1,1))
    plot(x)
    dev.off()
  }
  pdf(file=paste0(file, "_rates.pdf"))
  plotRates(x)
  dev.off()

#  if(!is.null(x$bs)) write.nexus(x$bs, file=paste0(file, "_bs.nex"),
#                                 digits=digits)
  sink(paste0(file, ".txt"))
  cat("phangorn", packageDescription("phangorn", fields = "Version"), "\n\n")
  print(summary(x$data))
  cat("\n\n")
  print(x)
  cat("\n\nThe following lines (re-)creates the pml object up to numerical inaccuracies:\n\n")
  call <- x$call
  call$data <- quote(align)
  call$tree <- quote(tree)
  cat("tree <- read.tree(\"", file, "_tree.nwk\")\n", sep="")
  type <- attr(x$data, "type")
  cat("align <- read.phyDat(\"", file,
      "_align.fasta\", format=\"fasta\", type=\"", type,"\")", sep="")
  cat( "\nfit <- ")
  print(call)
  if(save_rds){
    cat("\nAnd the following reproduces the exact pml object:\n\n")
    cat("fit <- readRDS(\"", file,".rds\")", sep = "")
  }
  cat("\n\nREFERENCES\n\n")
  cat("To cite phangorn please use:\n\n")
  print(citation("phangorn") [[1]], style = "text")
  sink()
  invisible(x)
}
