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
  if (x$inv > 0) cat("Proportion of invariant sites:", x$inv, "\n")
  if (x$k > 1) {
    cat("Model of rate heterogeneity: ")
    if(x$site.rate=="gamma") cat("Discrete gamma model\n")
    if(x$site.rate=="free_rate") cat("Free rate model\n")
    if(x$site.rate=="gamma_quadrature") cat("Discrete gamma model (quadrature) \n")
    if(x$site.rate=="gamma_unbiased") cat("Discrete gamma model (phangorn) \n")
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
  if(!isTRUE(all.equal(x$rate, 1))) cat("\nRate:", x$rate, "\n")
}


#' Export pml objects
#'
#' \code{write.pml} writes out the ML tree and the model parameters.
#'
#' @param x an object of class ancestral.
#' @param file a file name. File endings are added.
#' @param save_rds logical, if TRUE saves the pml object as a rds file,
#' otherwise the alignment is saved as a fasta file.
#' @param ... Further arguments passed to or from other methods.
#' @returns \code{write.pml}  returns the input x invisibly.
#' @seealso \code{\link{ancestral.pml}}, \code{\link{plotAnc}}
#' @examples
#' data(woodmouse)
#' fit <- pml_bb(woodmouse, "JC", rearrangement = "none")
#' write.pml(fit, "woodmouse")
#' unlink(c("woodmouse_pml.txt", "woodmouse_tree.nwk", "woodmouse.rds"))
#' @importFrom utils citation
#' @export
write.pml <- function(x, file="pml", save_rds=TRUE, ...){
  digits <- -1
  if (hasArg("digits")) digits <- list(...)$digits
  write.tree(x$tree, file=paste0(file, "_tree.nwk"))
  if(save_rds) saveRDS(x, file=paste0(file, ".rds"))
  else write.phyDat(x$data, file=paste0(file, "_align.fasta"), format="fasta")
  if(!is.null(x$bs)) write.nexus(x$bs, file=paste0(file, "_bs.nex"),
                                 digits=digits)
  sink(paste0(file, ".txt"))
  cat("phangorn", packageDescription("phangorn", fields = "Version"), "\n\n")
  print(x)
  cat("\n\n")
  cat("You can (re-)create the pml object using:\n\n")
  if(save_rds){
    cat("fit <- readRDS(\"", file,".rds\")", sep="")
  }
  else {
    call <- x$call
    call$data <- quote(align)
    call$tree <- quote(tree)
    cat("tree <- read.tree(\"", file, "_tree.nwk\")\n", sep="")
    cat("align <- read.phyDat(\"", file, "_align.fasta\", format=\"fasta\")",
        sep="")
    cat( "\nfit <- ")
    print(call)
  }
  cat("\n\nREFERENCES\n\n")
  cat("To cite phangorn please use:\n\n")
  print(citation("phangorn") [[1]], style="text")
  sink()
  invisible(x)
}
