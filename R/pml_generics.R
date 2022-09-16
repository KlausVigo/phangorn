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


#' @export
plot.pml <- function(x, type="phylogram", ...){
  type <- match.arg(type, c("phylogram","cladogram", "fan", "unrooted",
                            "radial", "tidy"))
  plot.phylo(x$tree, type=type, ...)
  if(is.rooted(x$tree) & (type %in% c("phylogram","cladogram"))) axisPhylo()
  else add.scale.bar()
}


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
    cat("Discrete gamma model\n")
    cat("Number of rate categories:", x$k, "\n")
    cat("Shape parameter:", x$shape, "\n")
  }
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
}
