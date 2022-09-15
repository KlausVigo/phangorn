#' Likelihood of a tree.
#'
#' \code{pml_bb} for pml black box infers a phylogentic tree infers a tree using
#' maximum likelihood (ML).
#'
#' \code{pml_bb} is a convenience function combining \code{pml} and
#' \code{optim.pml}. If no tree is supplied, the function will generate a
#' starting tree. If a modelTest object is supplied the model will be chosen
#' according to BIC.
#'
#' Currently very experimental and likely to change.
#' @param x An alignment of class (either class \code{phyDat}, \code{DNAbin} or
#' \code{AAbin}) or an object of class \code{modelTest}.
#' @param model A string providing model (e.g. "GTR+G(4)+I"). Not necessary if
#' a modelTest object is supplied.
#' @param method One of "unrooted", "ultrametric" or "tiplabeled".
#' @param rearrangement Type of tree tree rearrangements to perform, one of
#' "none", "NNI", "stochastic" or "ratchet"
#' @param start A starting tree can be supplied.
#' @param tip.dates A named vector of sampling times associated to the tips /
#' sequences.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{pml_bb} returns an object of class pml.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{optim.pml}}, \code{\link{modelTest}}, \code{\link{rtt}}
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(woodmouse)
#' tmp <- pml_bb(woodmouse)
#'
#' data(Laurasiatherian)
#' mt <- modelTest(Laurasiatherian)
#' fit <- pml_bb(mt)
#'}
#' @rdname pml_bb
#' @export
pml_bb <- function(x, model=NULL, rearrangement="stochastic",
         method="unrooted", start=NULL, tip.dates=NULL, ...){
  fit <- NULL
  type <- NULL
  method <- match.arg(method, c("unrooted", "ultrametric", "tipdated"))

  optRooted <- FALSE
  optRate <- FALSE
  if(method=="ultrametric" || method=="tipdated") optRooted <- TRUE
  if(method=="tipdated") optRate <- TRUE

  if(inherits(x, "AAbin") || inherits(x, "DNAbin")) x <- as.phyDat(x)
  if(inherits(x, "modelTest")){
    mt <- x
    fit <- as.pml(x)
    model <- x$Model[which.min(x[, "BIC"])]
  }
  if(inherits(x, "pml")){
    fit <- x
    if(is.null(model)) model <- guess_model(fit)
  }
  # eps=1e-7 10 * tau
  if(inherits(x, "phyDat")){
    if(is.null(start)) start <- candidate_tree(x, method=method,
                                               tip.dates = tip.dates, eps=1e-7)
    type <- attr(x, "type")
    if(is.null(model)){
      stop("Please supply a model!")
    } else {
      para <- split_model(x=model, type=type)
      fit <- pml(start, x, k=para$k, ASC=para$ASC)
    }
    if(method=="tipdated" && !is.null(attr(start, "rate")))
      fit <- update(fit, rate=attr(start, "rate"))
  }
  type <- attr(fit$data, "type")
  para <- split_model(model, type)
  if(type=="AA" && para$optFreq){
    fit <- optim.pml(fit, model=para$model, optGamma=para$optGamma,
                  optInv=para$optInv, optBf=TRUE, rearrangement=rearrangement,
                  optRate=optRate, optRooted=optRooted, ...)
  } else {
    fit <- optim.pml(fit, model=para$model, optGamma=para$optGamma,
                     optInv=para$optInv, rearrangement = rearrangement,
                     optRate=optRate, optRooted=optRooted, ...)
  }
  fit
}



##  check models
#' @export pml
split_model <- function(x="GTR + G(4) + I", type="DNA"){
  mods <- NULL
  if(type=="DNA") mods <- .dnamodels
  if(type=="AA") mods <- .aamodels
  if(type=="USER") mods <- .usermodels

  m <- strsplit(x, "\\+")[[1]]
  m <- trimws(m) # |> toupper()

  tmp <-  match(m, mods)
  if(all(is.na(tmp))) stop("Could not find model!")
  else pos <- tmp[!is.na(tmp)]
  if(length(pos)>1)  stop("Error, fould several models!")

  model <- mods[pos]
  m <- m[is.na(tmp)]

  optInv <- FALSE
  optGamma <- FALSE
  k <- 1L
  optFreq <- FALSE
  ASC <- FALSE
  if(length(m)>0){
    pos <- grep("G\\(", m)
    if(length(pos)==1){
      optGamma <- TRUE
      k_tmp <- sub("G\\(", "", m[pos])
      k_tmp <- sub("\\)", "", k_tmp)
      k <- as.integer(k_tmp)
      m <- m[-pos]
    }
  }
  if(length(m)>0){
    pos <- grep("I", m)
    if(length(pos)==1){
      optInv <- TRUE
      m <- m[-pos]
    }
  }
  if(length(m)>0){
    pos <- grep("F", m)
    if(length(pos)==1){
      optFreq <- TRUE
      m <- m[-pos]
    }
  }
  if(length(m)>0){
    pos <- grep("ASC", m)
    if(length(pos)==1){
      ASC <- TRUE
      m <- m[-pos]
    }
  }
  list(model=model, optFreq=optFreq, optInv=optInv, optGamma=optGamma, k=k,
       ASC=ASC)
}
