#' Likelihood of a tree.
#'
#' \code{pml_bb} for pml black box infers a phylogentic tree infers a tree using
#' maximum likelihood (ML).
#'
#' It combines the \code{pml}, \code{optim.pml}, \code{modelTest} and
#' \code{bootstrap.pml}. It first infers the best model, than optimizes this
#' model and last but not least performs bootstrap analysis.
#'
#' Currently very experimental and likly to change.
#' @param x An alignment of class (either class \code{phyDat}, \code{DNAbin} or
#' \code{AAbin}) or an object of class \code{modelTest}.
#' @param bs Number of bootstrap samples.
#' @param model A string providing model (e.g. "GTR+G(4)+I") otherwise
#' \code{\link{modelTest}} is called.
#' @param method Estimate "unrooted", "ultrametric" or "tiplabeled" tree.
#' @param rearrangement Type of tree tree rearrangements to perform, one of
#' "none", "NNI", "stochastic" or "ratchet"
## @param start A starting tree can be supplied.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{pml_bb} returns a list. It contains an object of class
#' \code{pml}, the best tree found so far and possibly objects of class
#' \code{modelTest} and the bootstrapped trees.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{optim.pml}}, \code{\link{bootstrap.pml}},
#' \code{\link{modelTest}},
#' @keywords cluster
#' @examples
#'
#' \dontrun{
#' data(woodmouse)
#' tmp <- pml_bb(woodmouse)
#'}
#' @rdname pml_bb
#' @export
pml_bb <- function(x, model=NULL, rearrangement="stochastic",
                   method=c("unrooted", "ultrametric", "tip.dated"),
                   bs=100, ...){ # start=NULL,
  mt <- NULL
  fit <- NULL
  type <- NULL

  if(inherits(x, "AAbin") || inherits(x, "DNAbin")) x <- as.phyDat(x)
  if(inherits(x, "modelTest")){
    mt <- x
    fit <- as.pml(x)
  }
  if(inherits(x, "phyDat")){
     if(is.null(model)){
       mt <- modelTest(x)
       fit <- as.pml(mt)
     }
     else {
       para <- split_model(x=model, type="DNA")
       if(is.null(start)) start <- candidate_tree(x, method=method)
       fit <- pml(start, x, k=para$k)
     }
  }
  if(is.null(model)) model <- guess_model(fit)
  type <- attr(fit$data, "type")
  model_terms <- split_model(model, type)

  #fit <- pml(fit, model=para$model, k=para$k)

  fit <- optim.pml(fit, model=model, optGamma=para$optGamma, optInv=para$optInv,
            optBf=para$optBf, rearrangement = rearrangement)

  tree <- fit$tree

  if(bs > 1){
    bs <- bootstrap.pml(fit, bs=bs, rearrangement="NNI", ...)
    tree <- plotBS(tree, bs, type="none")
  }
  list(fit=fit, modelTest=mt, bootstrap=bs, tree=tree)
}


##  check models
#' @export pml
split_model <- function(x="GTR + G(4) + I", type="DNA"){
  mods <- NULL
  if(type=="DNA") mods <- .dnamodels
  if(type=="AA") mods <- .aamodels
  #  if(type="USER") MK, MKv , SYM ER GTR

  m <- strsplit(x, "\\+")[[1]]
  m <- trimws(m) |> toupper()

  tmp <-  match(m, mods)
  if(all(is.na(tmp))) stop("Could not find model!")
  else pos <- tmp[!is.na(tmp)]
  if(length(pos)>1)  stop("Error, fould several models!")

  model <- mods[pos]
  m <- m[is.na(tmp)]

  optInv <- FALSE
  optGamma <- FALSE
  k=1L
  optFreq <- FALSE

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
  list(model=model, optFreq=optFreq, optInv=optInv, optGamma=optGamma, k=k)
}



rell <- function(x, B = 1000){
  weight <- as.integer(attr(x, "weight"))
  lw <- attr(x, "nr")
  X <- matrix(NA_integer_, lw, B)
  wvec <- rep( seq_len(lw), weight)
  for (i in 1:B) X[,i] <- tabulate(sample(wvec, replace = TRUE), nbins = lw)
  X
}




