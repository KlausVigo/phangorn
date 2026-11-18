#' @rdname pml
#' @aliases optim.pml
#' @export
optim.pml <- function(object, optNni = FALSE, optBf = FALSE, optQ = FALSE,
                      optInv = FALSE, optGamma = FALSE, optEdge = TRUE,
                      optRate = FALSE, optRooted = FALSE,
                      control = pml.control(), model = NULL,
                      rearrangement = ifelse(optNni, "NNI", "none"),
                      subs = NULL, ratchet.par = ratchet.control(), ...) {
  aLRT <- FALSE
  rearrangement <- match.arg(rearrangement,
                             c("none", "NNI", "ratchet", "stochastic", "multi2di"))
  optNni <- rearrangement !=  "none"
  perturbation <- rearrangement %in%
    c("ratchet", "stochastic", "multi2di")
  if(rearrangement == "ratchet") fbs <- vector("list", ratchet.par$minit)
  extras <- match.call(expand.dots = FALSE)$...
  pmla <- c("wMix", "llMix")
  wMix <- object$wMix
  llMix <- object$llMix
  ASC <- object$ASC
  site.rate <- object$site.rate
  optFreeRate <- FALSE
  if(site.rate=="free_rate" && optGamma){
    optFreeRate <- TRUE
    optGamma <- FALSE
  }
  data <- object$data
  type <- attr(data, "type")
  optModel <- FALSE
  if (is.null(model)) model <- object$model
  else{
    if (type == "AA") object <- update(object, model = model)
    optModel <- TRUE
  }
  if (is.null(llMix)) llMix <- 0
  if (!is.null(extras)) {
    names(extras) <- pmla[pmatch(names(extras), pmla)]
    existing <- match(pmla, names(extras))
    if (!is.na(existing[1]))
      wMix <- eval(extras[[existing[1]]], parent.frame())
    if (!is.na(existing[2]))
      llMix <- eval(extras[[existing[2]]], parent.frame())
  }
  tree <- object$tree
  call <- object$call

  addTaxa <- FALSE
  trace <- control$trace
  tau <- control$tau
  method <- "unrooted"
  is_ultrametric <- FALSE
  timetree <- FALSE
  stopifnot("rooted / unrooted tree needs at least 2 / 3 tips" =
              length(tree$tip.label) >= (2 + !optRooted))
  if (is.rooted(tree)) {
    if (!optRooted && optEdge) {
      tree <- unroot(tree)
      attr(tree, "order") <- NULL
      tree <- reorder(tree, "postorder")
      message("I unrooted the tree")
    }
    else{
      is_ultrametric <- is.ultrametric(tree, option=2)
      if(!is_ultrametric) {
        timetree <- TRUE
        method <- "tipdated"
        tip.dates <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
        tip.dates <- tip.dates - min(tip.dates)
      } else {
        method <- "ultrametric"
      }
    }
  }
  if (optNni) {
    if(!timetree){
      mapping <- map_duplicates(data)
      if (!is.null(mapping)) {
        orig.data <- data
        addTaxa <- TRUE
        tree <- drop.tip(tree, mapping[, 1])
        #
        if(method=="unrooted" && Ntip(tree)==2){
          tree <- midpoint(tree)
          optRooted <- TRUE
          method <- "ultrametric"
          message("I rooted the tree")
        }
        if(method=="unrooted") tree <- unroot(tree)
        attr(tree, "order") <- NULL
        tree <- reorder(tree, "postorder")
      }
    }
    if (!is.binary(tree)) tree <- multi2di(tree)
    optEdge <- TRUE
  }
  if (length(tree$tip.label) < (3 + !optRooted)) {
    optNni <- FALSE
    perturbation <- FALSE
  }
  stopifnot("rooted / unrooted tree needs at least 2 / 3 tips" =
              length(tree$tip.label) >= (2 + !optRooted))

  tree <- reorder(tree, "postorder")
  if (any(tree$edge.length < tau)) {
    tree <- minEdge(tree, tau, is_ultrametric)
    object <- update.pml(object, tree = tree)
  }
  if (optEdge && optRate && !timetree) {
    warning("You can't optimize edges and rates at the same time, only edges are optimized!", call. = FALSE)
    optRate <- FALSE
  }
  # 1x if  optEdge <- optEdge || optRooted || optNNI
  if (optRooted) {
    assert_phylo(tree, is_rooted = TRUE)
    optEdge <- TRUE
  }
  data <- subset(data, tree$tip.label)
  #  type <- attr(data, "type")
  #  if (type == "AA") {
  #    if(!is.null(model)) object <- update(object, model = model)
  #    model <- object$model
  #  }
  if (type == "CODON") {
    .syn <- synonymous_subs(code=attr(data, "code"))
    .sub <- tstv_subs(code=attr(data, "code"))
    if (is.null(model)) model <- "codon1"
    model <- match.arg(model, c("codon0", "codon1", "codon2", "codon3", "YN98"))
    dnds <- object$dnds
    tstv <- object$tstv
    if (!is.null(model)) {
      if (model == "codon0") optQ <- FALSE
      else  optQ <- TRUE
      if (model == "YN98") optBf <- TRUE
    }
    bf_choice <- object$frequencies
    freq_df <- df_freq_codon(bf_choice)
    if (bf_choice=="F3x4") {
      bf <- F3x4(data)
      bf_codon <- bf_by_codon(data)
      object <- update.pml(object, bf = bf)
    }
  }
  nr <- as.integer(attr(data, "nr"))
  nc <- as.integer(attr(data, "nc"))
  if (type == "DNA" && optModel) {
    object <- update(object, model = model)
    model <- object$model
    tmp <- subsChoice(model, has_gap_state(data))
    optBf <- tmp$optBf
    optQ <- tmp$optQ
    subs <- tmp$subs
  }
  if (type == "USER" && optModel) {
    tmp <- subsChoice_USER(model, nc)
    subs <- tmp$subs
    object <- update(object, model = model)
    model <- object$model
    optBf <- tmp$optBf
    optQ <- tmp$optQ
  }
  if (optInv && object$inv == 0) {
    inv <- 2/3 * const_site(data) / sum(attr(data, "weight"))
    object <- update.pml(object, inv = inv)
  }
  if (optBf && control$statefreq == "empirical" && type != "CODON") {
    bf <- baseFreq(data)
    object <- update.pml(object, bf = bf)
  }
  Q <- object$Q
  if (is.null(subs) && optQ) subs <- c(seq_len(length(Q) - 1), 0)
  bf <- object$bf
  eig <- object$eig
  inv <- object$inv
  k <- object$k
  if (k == 1 && optGamma) {
    optGamma <- FALSE
    message("only one rate class, ignored optGamma")
  }
  if(ASC && optInv){
    optInv <- FALSE
    message('cannot estimate invariant sites and Mkv model, ignored optInv')
  }
  shape <- object$shape
  w <- object$w
  g <- object$g
  ll0 <- object$logLik
  INV <- object$INV
  ll.0 <- object$ll.0
  rate <- object$rate
  ll <- ll0
  ll1 <- ll0
  opti <- TRUE
  RELL <- NULL
  if(ratchet.par$rell && perturbation){
    RELL <- init_rell(data, B=ratchet.par$bs)
    RELL <- update_rell(RELL, object$siteLik, tree)
  }
  nTips <- as.integer(length(tree$tip.label))
  on.exit({
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q,
                   levels = attr(data, "levels"), inv = inv, rate = rate,
                   g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix,
                   wMix = wMix, site = TRUE, ASC=ASC)
    if (addTaxa) {
      tree <- add.tips(tree, tips = mapping[, 1], where = mapping[, 2],
                       edge.length = rep(0, nrow(mapping)))
      data <- orig.data
      if (!is.null(RELL)){
        bs <- RELL$bs
        for(i in seq_along(bs)){
          bs[[i]] <- add.tips(bs[[i]], tips = mapping[, 1],
                              where = mapping[, 2], edge.length = rep(0, nrow(mapping)))
        }
        RELL$bs <- bs
      }
    }
    df <- ifelse(optRooted, tree$Nnode, length(tree$edge.length))
    if(site.rate == "free_rate" && k > 1) df <- df + 2 * k - 3L
    if(site.rate == "variable_gamma" && k > 1) df <- df + k - 1L
    dfQ <- ifelse(is.null(subs), length(unique(Q)) - 1, max(subs))
    df <- switch(type,
                 DNA = df + (k > 1) + (optInv | (inv > 0)) + length(unique(bf)) - 1 + dfQ,
                 AA = df + (k > 1) + (optInv | (inv > 0)) +
                   optBf * (length(unique(bf)) - 1),
                 CODON = df + (k > 1) + (optInv | (inv > 0)) + freq_df + (dnds != 1) +
                   (tstv != 1),
                 USER = df + (k > 1) + (optInv | (inv > 0)) + length(unique(bf)) - 1 + dfQ)
    names(bf) <- NULL
    object <- list(logLik = tmp$loglik, inv = inv, k = k, shape = shape,
                   Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik,
                   weight = attr(data, "weight"),
                   g = g, w = w, eig = eig, data = data, model = model, method=method,
                   INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll, call = call, df = df,
                   wMix = wMix, llMix = llMix, ASC=ASC, site.rate=site.rate)
    if (type == "CODON") {
      object$dnds <- dnds
      object$tstv <- tstv
      object$frequencies <- bf_choice
    }
    class(object) <- "pml"
    extras <- pairlist(bf = bf, Q = Q, inv = inv, shape = shape, k=k,
                       rate = rate, model=model, g=g, w=w)[c(optBf, optQ, optInv,
                                                             optGamma, optGamma, optRate, optModel,
                                                             optFreeRate, optFreeRate || site.rate=="variable_gamma")]

    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (!all(existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    object$call <- call
    if(!is.null(RELL)){
      bs <- RELL$bs
      class(bs) <- "multiPhylo"
      bs <- .compressTipLabel(bs)
      object$bs <- bs
      spl <- as.splits(bs)
      object$tree <- addConfidences(object$tree, spl)
    }
    if(rearrangement=="ratchet"){
      ind_null <- sapply(fbs, is.null)
      fbs <- fbs[!ind_null]
      class(fbs) <- "multiPhylo"
      object$abs <-  object$bs
      object$bs <- fbs
    }
    pml.free()
    return(object)
  })
  pml.init(data, k)

  if (optEdge) {
    res <- opt_Edge(tree, data, rooted = optRooted, eig = eig, w = w, g = g,
                    bf = bf, inv=inv, rate = rate, ll.0 = ll.0, INV = INV,
                    llMix = llMix, wMix=wMix, ASC=ASC,
                    control = pml.control(epsilon = 1e-07, maxit = 10, trace = trace,
                                          tau = tau))
    #    if (res[[2]] > ll) {
    ll <- res[[2]]
    tree <- res[[1]]
    #    }
  }
  rounds <- 1
  while (opti) {
    if (optRate) {
      res <- optimRate(tree, data, rate = rate, inv = inv,
                       INV = INV, Q = Q, bf = bf, eig = eig, k = k,
                       shape = shape, w = w, ll.0 = ll.0)
      if (trace > 0)
        cat("optimize rate: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="rate",
                  site.rate=site.rate)
    }
    if (optBf && control$statefreq=="estimated") {
      if (type=="CODON" && bf_choice=="F3x4") res <- optimF3x4(tree, data,
                                                               bf_codon = bf_codon, inv = inv, Q = Q, w = w, g = g, INV = INV,
                                                               rate = rate, k = k, llMix = llMix, ASC=ASC)
      else res <- optimBf(tree, data, bf = bf, inv = inv, Q = Q, w = w, g = g,
                          INV = INV, rate = rate, k = k, llMix = llMix, wMix=wMix, ASC=ASC)
      if (trace > 0)
        cat("optimize base frequencies: ", ll, "-->", max(res[[2]], ll), "\n")
      if (res[[2]] > ll) {
        bf <- res[[1]]
        eig <- edQt(Q = Q, bf = bf)
        if (inv > 0) ll.0 <- as.matrix(INV %*% (bf * inv))
        if (wMix > 0) ll.0 <- ll.0 + llMix
        ll <- res[[2]]
      }
    }
    if (optQ) {
      if (type == "CODON") {
        ab <- c(tstv, dnds)
        res <- switch(model,
                      YN98 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
                                        syn = .syn, bf = bf, w = w, g = g, inv = inv,
                                        INV = INV, ll.0 = ll.0, rate = rate, k = k,
                                        ab = log(ab), optK = TRUE, optW = TRUE),
                      codon1 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
                                          syn = .syn, bf = bf, w = w, g = g, inv = inv,
                                          INV = INV, ll.0 = ll.0, rate = rate, k = k,
                                          ab = log(ab), optK = TRUE, optW = TRUE),
                      codon2 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
                                          syn = .syn, bf = bf, w = w, g = g, inv = inv,
                                          INV = INV, ll.0 = ll.0, rate = rate, k = k,
                                          ab = log(ab), optK = FALSE, optW = TRUE),
                      codon3 = optimCodon(tree, data, Q = length(.sub), subs = .sub,
                                          syn = .syn, bf = bf, w = w, g = g, inv = inv,
                                          INV = INV, ll.0 = ll.0, rate = rate, k = k,
                                          ab = log(ab), optK = TRUE, optW = FALSE))
        tmp <- res[[5]]
        m <- length(tmp)
        dnds <- tmp[m]

        if (m > 1) tstv <- tmp[1]
      }
      else{
        res <- optimQ(tree, data, Q = Q, subs = subs, bf = bf, w = w, g = g,
                      inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k,
                      llMix = llMix, wMix=wMix, ASC=ASC)
      }
      Q <- res[[1]]
      eig <- edQt(Q = Q, bf = bf)
      if (trace > 0) cat("optimize rate matrix: ", ll, "-->", res[[2]], "\n")
      ll <- res[[2]]
    }
    ### start sitewise
    if (optInv && !optFreeRate) {
      res <- optimInv(tree, data, inv = inv, INV = INV, Q = Q,
                      bf = bf, eig = eig, k = k, shape = shape, w=w, rate = rate,
                      llMix = llMix, wMix=wMix)
      if (trace > 0)
        cat("optimize invariant sites: ", ll, "-->", max(res[[2]], ll), "\n")
      updateRates(res, ll, rate, shape, k, inv, wMix, update="inv",
                  site.rate=site.rate)
      ll.0 <- as.matrix(INV %*% (bf * inv))
      if (wMix > 0) ll.0 <- ll.0 + llMix
    }
    if (optGamma) {
      if(site.rate=="variable_gamma") {
        res <- optimGammaPhangorn(tree, data, w=w, shape = shape,                                  inv = inv, INV = INV,
                                  bf = bf, eig = eig, ll.0 = ll.0, rate = rate,
                                  llMix = llMix, wMix=wMix, ASC=ASC)

        shape <- res[[1]][k+1]
        w <- res[[1]][1:k]
        g <- discrete.gamma(shape, length(w), w = w)

        if (inv > 0){
          w <- (1 - inv) * w
          g <- g / (1 - inv)
        }
        if (trace > 0)
          cat("optimize shape parameter: ", ll, "-->", max(res[[2]], ll), "\n")
        ll <- res[[2]]
      } else {
        res <- optimGammaInv(tree, data, shape = shape, k = k, inv = inv,
                             optInv = optInv,
                             INV = INV, Q = Q, bf = bf, eig = eig, ll.0 = ll.0,
                             rate = rate, llMix = llMix, wMix=wMix, ASC=ASC,
                             site.rate=site.rate)
        if (trace > 0)
          cat("optimize shape parameter: ", ll, "-->", max(res[[2]], ll), "\n")

        inv <- res[[1]][2]
        updateRates(res, ll, rate, shape, k, inv, wMix, update="shape",
                    site.rate=site.rate)
        if(optInv){
          ll.0 <- as.matrix(INV %*% (bf * inv))
          if (wMix > 0) ll.0 <- ll.0 + llMix
        }
      }
    }
    if (optFreeRate) {
      tmp_ll <- ll
      res_test <- optimFreeRate(tree, data, g = g, k = k, w = w, inv = inv,
                                INV = INV, bf = bf, eig = eig, optInv=optInv,
                                ll.0 = ll.0, rate = rate, timetree=timetree)
      if(res_test$value > ll) {
        w <- res_test$w
        g <- res_test$r
        tree <- res_test$tree
        ll <- res_test$value
        inv <- res_test$inv
        ll.0 <- as.matrix(INV %*% (bf * inv))
        if (wMix > 0) ll.0 <- ll.0 + llMix
        if (trace > 0) cat("optimize free rate parameters: ", tmp_ll, "-->",
                           ll, "\n")
      }
    }
    ### end sitewise
    if (optEdge) {
      res <- opt_Edge(tree, data, rooted = optRooted, eig = eig, w = w, g = g,
                      bf = bf, inv=inv, rate = rate, ll.0 = ll.0,
                      llMix = llMix, wMix=wMix, ASC=ASC,
                      control = pml.control(epsilon = 1e-08, maxit = 10,
                                            trace = trace, tau = tau))
      ll <- res[[2]]
      tree <- res[[1]]
    }
    epsR <- 1e-8
    if (optNni) {
      res <- opt_nni(tree, data, rooted=optRooted, iter_max=5, trace=trace,
                     ll=ll, w = w, g = g, eig = eig, bf = bf, inv=inv,
                     ll.0 = ll.0, INV = INV, llMix = llMix, wMix=wMix, ASC=ASC,
                     RELL=NULL,
                     control = list(eps=1e-08, maxit=3, trace=trace-1, tau=tau))
      ll <- res$logLik
      tree <- res$tree
      swap <- res$swap
      rounds <- 1
      if (swap == 0) optNni <- FALSE
    }
    if ( perturbation && !optNni ) {
      maxR <- ratchet.par$iter
      maxit <- ratchet.par$maxit
      minit <- ratchet.par$minit
      kmax <- 1
      i <- 1
      if((rearrangement == "stochastic" || rearrangement == "ratchet") && optRooted){
        dm <- dist.ml(data, bf=bf, Q=Q, exclude = "pairwise")
      }
      for(i in seq_len(maxit)){
        if(rearrangement == "stochastic"){
          tree2 <- .stochastic_fun(tree, dm, tau, method, tip.dates,
                                   ratchet.par)
        } else if(rearrangement == "ratchet"){
          tree2 <- .ratchet_fun(tree, data, rooted=optRooted, w = w, g = g,
                                eig = eig, bf = bf, inv=inv,
                                rate=rate, ll.0 = ll.0, INV = INV, llMix = llMix,
                                wMix=wMix, ASC=ASC,
                                control=list(eps=1e-08, maxit=3, trace=0,
                                             tau=tau))
          fbs[[i]] <- tree2
        } else if(rearrangement == "multi2di"){
          tree2 <- di2multi(tree, tol=10*tau, tip2root=TRUE)
          if(any(tabulate(tree2$edge[,1])>2)){
            tree2 <- multi2di(tree2)
            if(!optRooted) tree2 <- unroot(tree2)
            tree2 <- minEdge(tree2, tau)
            tree2 <- reorder(tree2, "postorder")
          }
          #else return(list(tree, ll))
          else{
            i <- maxit
            break()
          }
        }
        res <- opt_Edge(tree2, data, rooted=optRooted, eig=eig, w=w, g=g, bf=bf,
                        inv=inv, rate=rate, ll.0=ll.0, llMix = llMix, wMix=wMix,
                        ASC=ASC,
                        control = pml.control(epsilon = 1e-08, maxit = 10,
                                              trace = 0, tau = tau))
        ll2 <- res[[2]]
        tree2 <- res[[1]]
        swap <- 1
        #        ll2 <- pml.fit(tree2, data, bf, shape = shape, k = k, Q = Q,
        #          levels = attr(data, "levels"), inv = inv, rate = rate,
        #          g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0,
        #          llMix = llMix, wMix = wMix, site = FALSE, Mkv=Mkv)
        res <- opt_nni(tree2, data, rooted=optRooted, iter_max=25, trace=trace-2,
                       ll=ll2, w = w, g = g, eig = eig, bf = bf, inv=inv,
                       rate=rate, ll.0 = ll.0, INV = INV, llMix = llMix,
                       wMix=wMix, ASC=ASC, RELL=RELL, ll_current=ll,
                       control=list(eps=1e-08, maxit=3, trace=trace-3, tau=tau))
        if (res$logLik > (ll + epsR)) {
          tree <- res$tree
          ll <- res$logLik
          kmax <- 1
        }
        else kmax <- kmax + 1
        if(!is.null(RELL)) RELL <- res$RELL
        if (trace > 0) cat("\rRatchet iteration ", i, ", best pscore so far:", ll)
        if ( (kmax >= maxR) && (i >= minit)) break()
      }
      optNni <- TRUE
      perturbation <- FALSE
      rounds <- 1

      #     cat("\rIteration: ", i, ". Best parsimony score so far: ", mp, sep="")
      if (trace >= 0)cat("\n")
    }
    if ( perturbation==FALSE &&
         ((abs((ll1 - ll) / ll)  < control$eps) || rounds > control$maxit)){
      #      if(aLRT){
      #        support <-SH_quartet(tree, data, rooted=optRooted, iter_max=1,
      #                    trace=trace,
      #                    ll=ll, w = w, g = g, eig = eig, bf = bf, inv=inv,
      #                    ll.0 = ll.0, INV = INV, llMix = llMix, wMix=wMix, ASC=ASC,
      #                    RELL=NULL,
      #                    control = list(eps=1e-08, maxit=3, trace=trace-1, tau=tau))
      #      }
      opti <- FALSE
    }

    rounds <- rounds + 1
    ll1 <- ll
  }
}
