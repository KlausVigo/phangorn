#
# Maximum likelihood estimation
#
discrete.gamma <- function (alpha, k) 
{
    if (k == 1) return(1)
    quants <- qgamma((1:(k - 1))/k, shape = alpha, rate = alpha)
    diff( c(0, pgamma(quants * alpha, alpha + 1),1)) * k
}


optimQ <- function (tree, data, Q=rep(1,6), subs=rep(1,length(Q)), trace = 0, ...) 
{
    m = length(Q)
    n = max(subs)
    ab = numeric(n)
#    ab = log(Q[match(1:n, subs)])    
    for(i in 1:n) ab[i]=log(Q[which(subs==i)[1]])
    fn = function(ab, tree, data, m, n, subs,...) {
        Q = numeric(m)
        for(i in 1:n)Q[subs==i] = ab[i]
        pml.fit(tree, data, Q = exp(Q),...)# Q^2, ...)
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = 10, control = list(fnscale = -1, 
        maxit = 25, trace = trace), tree = tree, data = data, m=m, n=n, subs=subs,...)
    Q = rep(1, m)
    for(i in 1:n) Q[subs==i] = exp(res[[1]][i])
    res[[1]] = Q
    res
}    

  
optimCodon <- function (tree, data, Q=rep(1,1830), subs=rep(1,length(Q)), syn = rep(0, length(Q)), trace = 0L, ab = c(0,0), optK=TRUE, optW=TRUE, ...) 
{
    m = length(Q)
    n = 1L # max(subs)

    fn = function(ab, tree, data, m, n, subs, syn, optK, optW, ...) {
        Q = numeric(m)
        Q[subs==1] = 0 # transversion
        if(optK) Q[subs==2] = ab[1] # transition
        else Q[subs==2] = 0
        if(optW) Q[syn==1] = Q[syn==1] + ab[2] # ab[n+1] dnds
        Q[syn<0] = -Inf
        pml.fit(tree, data, Q = exp(Q),...)# Q^2, ...)
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = Inf, control = list(fnscale = -1, 
        maxit = 25, trace = trace), tree = tree, data = data, m=m, n=n, 
        subs=subs, syn=syn, optK=optK, optW=optW, ...)
    ab = exp(res[[1]])
    Q[subs==1] = 1 # transversion
    if(optK) Q[subs==2] = ab[1] # transition
    else{ 
        Q[subs==2] = 1
        ab[1] = 1 
        }  
    if(optW) Q[syn==1] = Q[syn==1] * ab[2] # dnds
    else ab[2] = 1
    Q[syn<0] = 0
    res[[5]] = ab
    res[[1]] = Q
    res
} 


subsChoice <- function(type=c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR")){
    type = match.arg(type)
    switch(type,
         JC = list(optQ=FALSE, optBf=FALSE,   subs=c(0, 0, 0, 0, 0, 0)),
         F81 = list(optQ=FALSE, optBf=TRUE,   subs=c(0, 0, 0, 0, 0, 0)),
         K80 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 0, 0, 1, 0)),
         HKY = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 1, 0)),
         TrNe = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 0, 0, 2, 0)),
         TrN = list(optQ=TRUE, optBf=TRUE,    subs=c(0, 1, 0, 0, 2, 0)),
         TPM1 = list(optQ=TRUE, optBf=FALSE,  subs=c(0, 1, 2, 2, 1, 0)),
         K81 = list(optQ=TRUE, optBf=FALSE,   subs=c(0, 1, 2, 2, 1, 0)),
         TPM1u = list(optQ=TRUE, optBf=TRUE,  subs=c(0, 1, 2, 2, 1, 0)),
         TPM2 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM2u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 1, 0, 2, 0)),
         TPM3 = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 0, 1, 2, 0)),
         TPM3u = list(optQ=TRUE, optBf=TRUE,  subs=c(1, 2, 0, 1, 2, 0)),
         TIM1e = list(optQ=TRUE, optBf=FALSE, subs=c(0, 1, 2, 2, 3, 0)),
         TIM1 = list(optQ=TRUE, optBf=TRUE,   subs=c(0, 1, 2, 2, 3, 0)),
         TIM2e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 1, 0, 3, 0)),
         TIM2 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 1, 0, 3, 0)),
         TIM3e = list(optQ=TRUE, optBf=FALSE, subs=c(1, 2, 0, 1, 3, 0)),
         TIM3 = list(optQ=TRUE, optBf=TRUE,   subs=c(1, 2, 0, 1, 3, 0)),
         TVMe = list(optQ=TRUE, optBf=FALSE,  subs=c(1, 2, 3, 4, 2, 0)),
         TVM = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 2, 0)),
         SYM = list(optQ=TRUE, optBf=FALSE,   subs=c(1, 2, 3, 4, 5, 0)),
         GTR = list(optQ=TRUE, optBf=TRUE,    subs=c(1, 2, 3, 4, 5, 0))
         )
}


modelTest <- function (object, tree = NULL, model = c("JC", "F81", "K80", 
    "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, k = 4, control = pml.control(epsilon = 1e-08, 
    maxit = 10, trace = 1), multicore = FALSE) 
{    
    if (class(object) == "phyDat") 
        data = object
    if (class(object) == "pml") {
        data = object$data
        if (is.null(tree)) 
            tree = object$tree
    }
    
    if(attr(data, "type")=="DNA") type = c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", 
            "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", 
            "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", 
            "SYM", "GTR")
    if(attr(data, "type")=="AA") type = .aamodels   
    model = match.arg(model, type, TRUE)
    
    env = new.env()
    assign("data", data, envir=env)
    
    if (is.null(tree)) 
        tree = NJ(dist.hamming(data))
    else{
        tree <- nnls.phylo(tree, dist.ml(data)) 
        # may need something faster for trees > 500 taxa  
    }
    trace <- control$trace
    control$trace = trace - 1
    fit = pml(tree, data)
    fit = optim.pml(fit, control = control)
    l = length(model)
    n = 1L + sum(I + G + (G & I))
    nseq = sum(attr(data, "weight"))
    fitPar = function(model, fit, G, I, k) {
        m = 1
        res = matrix(NA, n, 6)
        res = as.data.frame(res)
        colnames(res) = c("Model", "df", "logLik", "AIC", "AICc", "BIC")
        data.frame(c("Model", "df", "logLik", "AIC", "AICc", "BIC"))
        calls = vector("list", n)
        trees = vector("list", n)
        fittmp = optim.pml(fit, model = model, control = control)
        res[m, 1] = model
        res[m, 2] = fittmp$df
        res[m, 3] = fittmp$logLik
        res[m, 4] = AIC(fittmp)
        res[m, 5] = AICc(fittmp)
        res[m, 6] = AIC(fittmp, k = log(nseq))
        calls[[m]] = fittmp$call
        
        trees[[m]] = fittmp$tree
        m = m + 1
        if (I) {
            if(trace>0)print(paste(model, "+I", sep = ""))
            fitI = optim.pml(fittmp, model = model, optInv = TRUE, 
                             control = control)
            res[m, 1] = paste(model, "+I", sep = "")
            res[m, 2] = fitI$df
            res[m, 3] = fitI$logLik
            res[m, 4] = AIC(fitI)
            res[m, 5] = AICc(fitI)
            res[m, 6] = AIC(fitI, k = log(nseq))
            calls[[m]] = fitI$call
            trees[[m]] = fitI$tree
            m = m + 1
        }
        if (G) {
            if(trace>0)print(paste(model, "+G", sep = ""))
            fitG = update(fittmp, k = k)
            fitG = optim.pml(fitG, model = model, optGamma = TRUE, 
                             control = control)
            res[m, 1] = paste(model, "+G", sep = "")
            res[m, 2] = fitG$df
            res[m, 3] = fitG$logLik
            res[m, 4] = AIC(fitG)
            res[m, 5] = AICc(fitG)
            res[m, 6] = AIC(fitG, k = log(nseq))
            calls[[m]] = fitG$call
            trees[[m]] = fitG$tree
            m = m + 1
        }
        if (G & I) {
            if(trace>0)print(paste(model, "+G+I", sep = ""))
            fitGI = optim.pml(fitG, model = model, optGamma = TRUE, 
                              optInv = TRUE, control = control)
            res[m, 1] = paste(model, "+G+I", sep = "")
            res[m, 2] = fitGI$df
            res[m, 3] = fitGI$logLik
            res[m, 4] = AIC(fitGI)
            res[m, 5] = AICc(fitGI)
            res[m, 6] = AIC(fitGI, k = log(nseq))
            calls[[m]] = fitGI$call
            trees[[m]] = fitGI$tree
            m = m + 1
        }
        list(res, trees, calls)
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {
        # !require(parallel) ||         
        if (.Platform$GUI != "X11") {
            warning("package 'parallel' not found or GUI is used, \n      analysis is performed in serial")
        }
        else {
            RES <- mclapply(model, fitPar, fit, G, I, k)
            eval.success <- TRUE
        }
    }
    if (!eval.success) 
        res <- RES <- lapply(model, fitPar, fit, G, I, k)
    
    RESULT = matrix(NA, n * l, 6)
    RESULT = as.data.frame(RESULT)
    colnames(RESULT) = c("Model", "df", "logLik", "AIC", "AICc", "BIC")
    for (i in 1:l) RESULT[((i - 1) * n + 1):(n * i), ] = RES[[i]][[1]]
    for(i in 1:l){
        for(j in 1:n){
            mo = RES[[i]][[1]][j,1]
            tname = paste("tree_", mo, sep = "")
            tmpmod = RES[[i]][[3]][[j]]
            tmpmod["tree"] = call(tname)
            if(!is.null(tmpmod[["k"]]))tmpmod["k"] = k
            if(attr(data, "type")=="AA") tmpmod["model"] = RES[[i]][[1]][1,1]          
            assign(tname, RES[[i]][[2]][[j]], envir=env)
            assign(mo, tmpmod, envir=env) 
        }
    }
    attr(RESULT, "env") = env 
    RESULT
}

 
optimGamma = function(tree, data, shape=1, k=4,...){
    fn = function(shape, tree, data, k,...)pml.fit(tree, data, shape=shape, k=k,...)
    res = optimize(f=fn, interval = c(0.1, 500), lower = 0.1, upper = 500, maximum = TRUE,
        tol = .01, tree=tree, data=data, k=k,...)
    res
    }
    
 
optimInv = function(tree, data, inv=0.01, INV=NULL, ll.0=NULL,...){
    fn = function(inv, tree, data,...)pml.fit(tree, data, inv=inv, INV=INV, ll.0=NULL,...)
    res = optimize(f=fn, interval = c(0,1), lower = 0, upper = 1, maximum = TRUE,
         tol = .0001, tree=tree, data=data,...)
    res
    }
  

# changed to c(-10,10) from c(-5,5)
optimRate <- function(tree, data, rate=1, ...){
    fn <- function(rate, tree, data, ...) pml.fit(tree, data, rate=exp(rate), ...)
    res <- optimize(f = fn, interval = c(-10, 10), tree = tree, data = data, ..., maximum = TRUE)
    res[[1]] <- exp(res[[1]])
    res
}
    

optimBf = function(tree, data, bf=c(.25,.25,.25,.25), trace=0,...){
    l=length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, tree, data,...){
        bf = exp(c(lbf,0))
        bf = bf/sum(bf)
        pml.fit(tree, data, bf=bf, ...)
        }
    res = optim(par=lbf, fn=fn, gr=NULL, method="Nelder-Mead", control=list(fnscale=-1, maxit=500, trace=trace),tree=tree, data=data,...)
    bf = exp(c(res[[1]],0))
    bf = bf/sum(bf)
    result = list(bf=bf, loglik = res[[2]])
    result
    }


optimW = function(fit,...){
    w = fit$w
    g = fit$g
    siteLik = fit$siteLik
    k = length(w)
    l = dim(siteLik[[1]])[1]
    x=matrix(0,l,k)
    for(i in 1:k)x[,i] = rowSums(siteLik[[i]])
    weight = fit$weight
    nenner = 1/w[k]
    eta = log(w * nenner)
    eta = eta[-k]
    fn = function(eta,x,g,weight){
        eta = c(eta,0)
        p = exp(eta)/sum(exp(eta))
        res = x%*%p
        res = sum(weight*log(res))  * (1 + abs(sum(p*g) - 1))
        res
    }  
    res = optim(eta, fn = fn, method = "Nelder-Mead", control=list(fnscale=-1, reltol = 1e-12),gr=NULL, x=x,g=g, weight=weight)
    p = exp(c(res$par,0))
    p = p/sum(p)
    result = list(par = p, value = res$value)
    result    
}


#predict.pml <- function(object, newdata,...) sum(object$site * newdata)


logLik.pml <- function(object,...){
    res <- object$logLik
    attr(res,"df") <- object$df
    class(res) <- "logLik"
    res
}


AICc <- function (object, ...) 
    UseMethod("AICc")


AICc.pml <- function(object, ...){
    n = sum(object$weight)
    k = object$df
#    if(k>=(n-1))return(NULL)    
    res = AIC(object)
    res +   (2*k*(k+1))/(n-k-1)    
}


anova.pml <- function (object, ...) 
{
    X <- c(list(object), list(...))
    df <- sapply(X, "[[", "df")
    ll <- sapply(X, "[[", "logLik")
    dev <- c(NA, 2 * diff(ll)) 
    ddf <- c(NA, diff(df))
    table <- data.frame(ll, df, ddf, dev, pchisq(dev, ddf, lower.tail = FALSE))
    dimnames(table) <- list(1:length(X), c("Log lik.", "Df", 
        "Df change", "Diff log lik.", "Pr(>|Chi|)"))
    structure(table, heading = "Likelihood Ratio Test Table", 
        class = c("anova", "data.frame"))
}
    
    
#vcov.pml <- function(object, obs=FALSE,...){
#    if(obs) FI = score4(object)[[2]]
#    else FI = score(object,FALSE)[[2]]
#    l = dim(FI)[1]
#    res = try(solve(FI))
#    if(class(res) == "try-error"){
#        cat("Covariance is ill-conditioned !! \n")
#        res = solve(FI + diag(l)* 1e-8)
#        }
#    res
#}
                             
vcov.pml <- function(object, ...){
    FI = score(object,FALSE)[[2]]
    l = dim(FI)[1]
    res = try(solve(FI))
    if(class(res) == "try-error"){
        cat("Covariance is ill-conditioned !! \n")
        res = solve(FI + diag(l)* 1e-8)
        }
    res
}


getd2P <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getd2PM",eig,as.integer(n),as.double(el),as.double(g))
    attr(res,"dim") <- c(length(g),length(el))
    res
}


getdP <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getdPM",eig,as.integer(n),as.double(el),as.double(g))
    attr(res,"dim") <- c(length(g),length(el))
    res
}


# version without transformation (used for vcov)
getdP2 <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getdPM2",eig,as.integer(n),as.double(el),as.double(g))
    attr(res,"dim") <- c(length(g),length(el))
    res
}


# version without transformation 
getd2P2 <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)    
    res <- .Call("getd2PM2",eig,as.integer(n),as.double(el),as.double(g))
    attr(res,"dim") <- c(length(g),length(el))
    res
}


getP <- function(el, eig=edQt(), g=1.0){
    n <- length(eig$values)
    res <- .Call("getPM", eig, as.integer(n), as.double(el), as.double(g))
    attr(res, "dim") <- c(length(g), length(el)) 
    res
}


lli <- function (data, tree=NULL, ...) 
{
    contrast <- attr(data, "contrast")
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")
    nco <- as.integer(dim(contrast)[1])
    if(!is.null(tree)) data <- subset(data, tree$tip.label)
    .Call("invSites", data, as.integer(nr), as.integer(nc), contrast, as.integer(nco))    
}


edQt <- function (Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25)) 
{
    l = length(bf)
    res = matrix(0, l, l)
    res[lower.tri(res)] = Q
    res = res + t(res)
    res = res * bf
    res2 = res * rep(bf, each = l)    
    diag(res) = -colSums(res)
    res = res/sum(res2)
    e = eigen(res, FALSE)
    e$inv = solve.default(e$vec)
    e
}


edQ <- function(Q=c(1,1,1,1,1,1), bf=c(0.25,.25,.25,.25)){
    l=length(bf)
    res = matrix(0, l, l)
    res[lower.tri(res)] = Q
    res = res+t(res)
    res = res * rep(bf,each=l)
    diag(res) = -rowSums(res)
    res2 = res * rep(bf,l)
    diag(res2)=0 
    res = res/sum(res2)
    e = eigen(res, FALSE)
    e$inv = solve.default(e$vec)
    e
}

edQ2 <- function(Q){
    res = Q
    l=dim(Q)[1]
    diag(res) = 0
    diag(res) = -rowSums(res)
    e = eigen(res, FALSE)
    e$inv = solve.default(e$vec)
    e
}


pml.free <- function(){
    .C("ll_free")
#    rm(.INV, .iind, envir = parent.frame())
}


pml.init <- function(data, k=1L){
    nTips <- length(data)
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")    
    .C("ll_init", as.integer(nr), as.integer(nTips), as.integer(nc), as.integer(k))
    INV <- lli(data) #, tree
#    .iind <<- which((INV %*% rep(1, nc)) > 0)
#    .INV <<-  Matrix(INV, sparse=TRUE)
    assign(".iind", which((INV %*% rep(1, nc)) > 0), envir=parent.frame())
    assign(".INV", Matrix(INV, sparse=TRUE), envir=parent.frame())
} 



pml.free2 <- function(){.C("ll_free2")}

pml.init2 <- function(data, k=1L){
    nTips <- length(data)
    nr <- attr(data, "nr")
    nc <- attr(data, "nc")    
    weight <- attr(data, "weight")
    .C("ll_init2", as.integer(unlist(data, use.names=FALSE)), as.double(weight), as.integer(nr), as.integer(nTips), as.integer(nc), as.integer(k))
} 


fn.quartet <- function(old.el, eig, bf, dat,  g=1, w=1, weight, ll.0) {
    l= length(dat[,1]) 
    ll = ll.0
    res = vector("list", 2*l)
    tmp1 = NULL
    tmp2 = NULL
    attr(res,"dim") = c(l,2)
    for(j in 1:l){
            P = getP(old.el, eig, g[j])
            tmp1 = (dat[[j,1]] %*% P[[1]]) *(dat[[j,2]] %*% P[[2]])
            tmp2 = (dat[[j,3]] %*% P[[3]]) * (dat[[j,4]] %*% P[[4]])
            res[[j,1]] = tmp1 * (tmp2 %*% P[[5]])
            res[[j,2]] = tmp2
            ll = ll +  res[[j,1]] %*% (w[j]*bf)
        } 
    l0 = sum(weight * log(ll))
    list(ll=l0,res=res)
}


fn.quartet2 <- function (old.el, eig, bf, dat1, dat2, dat3, dat4, g = 1, w = 1, 
    weight, ll.0, contrast, ext) 
{
    l = length(w)
    ll = ll.0
    res = vector("list", 2 * l)
    tmp1 = NULL
    tmp2 = NULL
    attr(res, "dim") = c(l, 2)
    for (j in 1:l) {
        P = getP(old.el, eig, g[j])
        if (ext[1] == FALSE && ext[2] == FALSE) 
            tmp1 = (dat1[[j]] %*% P[[1]]) * (dat2[[j]] %*% P[[2]])
        if (ext[1] == FALSE && ext[2] == TRUE) 
            tmp1 = (dat1[[j]] %*% P[[1]]) * (contrast %*% P[[2]])[dat2, ]
        if (ext[1] == TRUE && ext[2] == FALSE) 
            tmp1 = (contrast %*% P[[1]])[dat1, ] * (dat2[[j]] %*% P[[2]])
        if (ext[1] == TRUE && ext[2] == TRUE) 
            tmp1 = (contrast %*% P[[1]])[dat1, ] * (contrast %*% P[[2]])[dat2, ]
        if (ext[3] == FALSE && ext[4] == FALSE) 
            tmp2 = (dat3[[j]] %*% P[[3]]) * (dat4[[j]] %*% P[[4]])
        if (ext[3] == FALSE && ext[4] == TRUE) 
            tmp2 = (dat3[[j]] %*% P[[3]]) * (contrast %*% P[[4]])[dat4, ]
        if (ext[3] == TRUE && ext[4] == FALSE) 
            tmp2 = (contrast %*% P[[3]])[dat3, ] * (dat4[[j]] %*% P[[4]])
        if (ext[3] == TRUE && ext[4] == TRUE) 
            tmp2 = (contrast %*% P[[3]])[dat3, ] * (contrast %*% P[[4]])[dat4, ]
        res[[j, 1]] = tmp1 * (tmp2 %*% P[[5]])
        res[[j, 2]] = tmp2
        ll = ll + res[[j, 1]] %*% (w[j] * bf)
    }
    l0 = sum(weight * log(ll))
    list(ll = l0, res = res)
}


optim.quartet2 <- function (old.el, eig, bf, dat1, dat2, dat3, dat4, g = 1, w = 1, 
    weight, ll.0 = weight * 0, control = list(eps = 1e-08, maxit = 5, 
        trace = 0), llcomp = -Inf, evi, contrast, contrast2, 
    ext = c(FALSE, FALSE, FALSE, FALSE)) 
{
    eps = 1
    iter = 0
    while (eps > control$eps && iter < control$maxit) {
        tmp <- fn.quartet2(old.el = old.el, eig = eig, bf = bf, 
            dat1 = dat1, dat2 = dat2, dat3 = dat3, dat4 = dat4, 
            g = g, w = w, weight = weight, ll.0 = ll.0, contrast=contrast, ext = ext)
        old.ll = tmp$ll
   
        el1 <- fs3(old.el[1], eig, tmp$res[, 1], dat1, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[1], getA=TRUE, getB=FALSE)
        el2 <- fs3(old.el[2], eig, el1[[2]], dat2, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[2], getA=TRUE, getB=FALSE)
        el5 <- fs3(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = 0L, getA=FALSE, getB=TRUE)
        el3 <- fs3(old.el[3], eig, el5[[3]], dat3, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[3], getA=TRUE, getB=FALSE)
        el4 <- fs3(old.el[4], eig, el3[[2]], dat4, weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, contrast=contrast, contrast2=contrast2, evi=evi, ext = ext[4], getA=FALSE, getB=FALSE)
        old.el[1] = el1[[1]]
        old.el[2] = el2[[1]]
        old.el[3] = el3[[1]]
        old.el[4] = el4[[1]]
        old.el[5] = el5[[1]]
        iter = iter + 1
        ll = el4[[4]]
        eps = (old.ll - ll)/ll
        if (ll < llcomp) 
            return(list(old.el, ll))
        old.ll = ll
    }
    list(old.el, ll)
}


pml.nni <- function (tree, data, w, g, eig, bf, ll.0, ll, ...) 
{        
    k = length(w)
    INDEX <-  indexNNI(tree)
    rootEdges <- attr(INDEX,"root")
    .dat <- NULL
    data = getCols(data, tree$tip)

    parent = tree$edge[,1]
    child = tree$edge[,2]
    weight = attr(data, "weight")
    datp = rnodes(tree, data, w, g, eig, bf)    
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)

    nTips = length(tree$tip.label)
    evector <- numeric(max(parent)) 
    evector[child] <- tree$edge.length
    m <- dim(INDEX)[1]
    loglik = numeric(2*m)
    edgeMatrix <- matrix(0, 2*m, 5)
    l = length(datp[, 1])
    for(i in 1:m){
        ei = INDEX[i,]
        el0 = evector[INDEX[i,]]

        ext = ei[1:4] < nTips+1L
        if (!(ei[5] %in% rootEdges)) dat1 = datp[, ei[1], drop = FALSE]
        else{ if(ext[1]) dat1 = data[[ ei[1] ]]
             else dat1 = .dat[, ei[1], drop=FALSE]
        } 
        if(ext[2]) dat2 = data[[ ei[2] ]]
             else dat2 = .dat[, ei[2], drop=FALSE] 
        if(ext[3]) dat3 = data[[ ei[3] ]]
             else dat3 = .dat[, ei[3], drop=FALSE]
        if(ext[4]) dat4 = data[[ ei[4] ]]
             else dat4 = .dat[, ei[4], drop=FALSE]

        new1 <- optim.quartet2(el0[c(1, 3, 2, 4, 5)], eig, bf, 
            dat1, dat3, dat2, dat4, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 3, 2, 4)])
        new2 <- optim.quartet2(el0[c(1, 4, 3, 2, 5)], eig, bf,  
            dat1, dat4, dat3, dat2, g, w, weight, ll.0, llcomp=ll, evi=evi, contrast=contrast, contrast2=contrast2, ext=ext[c(1, 4, 3, 2)])


        loglik[(2*i)-1]=new1[[2]]
        loglik[(2*i)]=new2[[2]] 
        edgeMatrix[(2*i)-1,]=new1[[1]]
        edgeMatrix[(2*i),]=new2[[1]]           
    }
    swap <- 0
    eps0 <- 1e-6
    candidates <- loglik > ll + eps0

    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label))
 
#    on.exit(.C("ll_free"))
#    .C("ll_init", nr, nTips, nc, as.integer(k))

    while(any(candidates)){     
        ind = which.max(loglik)
        loglik[ind]=-Inf
        if( ind %% 2 ) swap.edge = c(2,3)
        else swap.edge = c(2,4)
        tree2 <- changeEdge(tree, INDEX[(ind+1)%/%2,swap.edge], INDEX[(ind+1)%/%2,], edgeMatrix[ind,])
 
        test <- pml.fit(tree2, data, bf = bf, k=k, g=g, w=w, eig=eig, ll.0=ll.0, ...) 
        if(test <= ll + eps0) candidates[ind] = FALSE
        if(test > ll + eps0) {
            ll = test 
            swap=swap+1
            tree <- tree2
            indi <- which(rep(colSums(apply(INDEX,1,match,INDEX[(ind+1)%/%2,],nomatch=0))>0,each=2))
            candidates[indi] <- FALSE
            loglik[indi] <- -Inf
        }
    } 
    list(tree=tree, ll=ll, swap=swap)     
}


rnodes <- function (tree, data, w, g, eig, bf) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    data = getCols(data, tree$tip) 
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    l = length(w)        
    dat = vector(mode = "list", length = m*l)
    dim(dat) <- c(l,m)
    tmp = length(data)
#    for(i in 1:length(w))dat[i,1:tmp]=new2old.phyDat(data) #
#    dat[1,1:tmp] <- data  vielleicht gebraucht
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    for(i in 1:l)dat[i,(q + 1):m] <- .Call("LogLik2", data, P[i,], nr, nc, node, edge, nTips, mNodes, contrast, nco)
    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    datp = vector("list", m)   
    dat2 = vector("list", m * l)
    dim(dat2) <- c(l,m)
    for(i in 1:l){     
      datp[(nTips + 1)] = dat[i,(nTips + 1)]
      for (j in (m - 1):1) {
          if (child[j] > nTips){
             tmp2 = (datp[[parent[j]]]/(dat[[i,child[j]]] %*% P[[i,j]]))
             datp[[child[j]]] = (tmp2 %*% P[[i,j]]) * dat[[i,child[j]]]  
             dat2[[i, child[j]]] = tmp2
             }
       }
    }
    assign(".dat", dat, envir = parent.frame(n = 1))
    dat2
}


score <- function (fit, transform=TRUE) 
{
    tree = fit$tree
    child <- tree$edge[, 2]
    l = length(child)
    sc = numeric(l)
    weight = as.numeric(fit$weight)
    f <- drop(exp(fit$site))
    dl = dl(fit, transform)
    dl = dl/f
    sc = colSums(weight * dl)
    F = crossprod(dl*weight,dl) 
    names(sc) = child
    dimnames(F) = list(child, child) 
    result = list(sc = sc, F = F)
    result
}


# wird noch in partition models verwendet
optim.quartet <- function (old.el, eig, bf, dat, g = 1, w = 1, weight, ll.0 = weight * 
    0, control = list(eps = 1e-08, maxit = 5, trace = 0), llcomp=-Inf) 
{
    eps = 1
    iter = 0
    evi = (t(eig[[3]]) * bf)
    while (eps > control$eps && iter < control$maxit) {
        tmp <- fn.quartet(old.el = old.el, eig = eig, bf = bf, dat = dat, 
            g = g, w = w, weight = weight, ll.0 = ll.0)
        old.ll = tmp$ll 
        el1 <- fs(old.el[1], eig, tmp$res[, 1], dat[, 1], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el2 <- fs(old.el[2], eig, el1[[2]], dat[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el5 <- fs(old.el[5], eig, el2[[2]], tmp$res[, 2], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=FALSE, getB=TRUE)
        el3 <- fs(old.el[3], eig, el5[[3]], dat[, 3], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=TRUE, getB=FALSE)
        el4 <- fs(old.el[4], eig, el3[[2]], dat[, 4], weight, 
            g = g, w = w, bf = bf, ll.0 = ll.0, evi, getA=FALSE, getB=FALSE)
        old.el[1] = el1[[1]]
        old.el[2] = el2[[1]]
        old.el[3] = el3[[1]]
        old.el[4] = el4[[1]]
        old.el[5] = el5[[1]]
        iter = iter + 1
        ll = el4[[4]]
        eps = (old.ll - ll) / ll
        if(ll<llcomp)return(list(old.el, ll))  
        old.ll = ll
    }
    list(old.el, ll)
}


plot.pml<-function(x,...)plot.phylo(x$tree,...)


phangornParseFormula <- function(model){

    parseSide <- function(model) {
        model.vars <- list()
        while (length(model) == 3 && model[[1]] == as.name("+")) {
            model.vars <- c(model.vars, model[[3]])
            model <- model[[2]]
        }
        unlist(rev(c(model.vars, model)))

    } 

    if (!inherits(model, "formula")) 
        stop("model must be a formula object")
    l <- length(model)
    varsLHS <- NULL       
    if(l==3){        
        modelLHS <- model[[2]]
        modelRHS <- model[[3]]
        varsRHS <- parseSide(modelRHS)
        varsRHS <- unlist(lapply(varsRHS,as.character))
        varsLHS <- parseSide(modelLHS)
        varsLHS <- unlist(lapply(varsLHS,as.character))
    }
    if(l==2){
       modelRHS <- model[[2]]
       varsRHS <- parseSide(modelRHS)
       varsRHS <- unlist(lapply(varsRHS,as.character))
    }
    list(left=varsLHS, right=varsRHS)
}


pml.control <- function (epsilon = 1e-08, maxit = 10, trace = 1) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}


optim.pml <- function (object, optNni = FALSE, optBf = FALSE, optQ = FALSE, 
    optInv = FALSE, optGamma = FALSE, optEdge = TRUE, optRate = FALSE, optRooted=FALSE, 
    control = pml.control(epsilon = 1e-8, maxit = 10, trace = 1L), 
    model = NULL, subs = NULL, ...) 
{
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix")
    wMix <- object$wMix
    llMix <- object$llMix
    if(is.null(llMix)) llMix=0
    if (!is.null(extras)) {
        names(extras) <- pmla[pmatch(names(extras), pmla)]
        existing <- match(pmla, names(extras))
        if (!is.na(existing[1])) 
            wMix <- eval(extras[[existing[1]]], parent.frame())
        if (!is.na(existing[2])) 
            llMix <- eval(extras[[existing[2]]], parent.frame())
    }
    tree = object$tree
    call = object$call
    if(optNni) {
        if(!is.binary.tree(tree)) 
            tree = multi2di(tree)
        optEdge = TRUE     
    }
    if(is.rooted(tree)) {
        if(optRooted==FALSE && optEdge==TRUE){
            tree = unroot(tree)
            attr(tree, "order") <- NULL
            tree = reorder(tree, "postorder")
            warning("I unrooted the tree", call. = FALSE)
        }    
    }
    if(is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    if(any(tree$edge.length < 1e-08)) {
        tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
# save to change to new update.pml       
        object <- update.pml(object, tree = tree)
    }
    if(optEdge & optRate) {
        warning("You can't optimise edges and rates at the same time, only edges are optimised!", call. = FALSE)
        optRate = FALSE
    }
    if(optRooted){
        optEdge = FALSE
        if(!is.rooted(tree)) stop("Tree must be rooted!")
        if(!is.ultrametric(tree)) stop("Tree must be ultrametric!")
	}
    trace <- control$trace
    
    data = object$data
    data = subset(data, tree$tip.label) 

    type <- attr(data, "type")
    if (type == "AA" & !is.null(model)){
        object = update(object, model=model)  
    }     
    if (type == "CODON") {
        dnds <- object$dnds 
        tstv <- object$tstv
        if(!is.null(model)){
            if(model == "codon0") optQ = FALSE
            else  optQ = TRUE
        }
    }       
    Q = object$Q
    if(is.null(subs)) subs = c(1:(length(Q) - 1), 0)
    bf = object$bf
    eig = object$eig
    inv = object$inv
    k = object$k
    if(k==1 & optGamma){
        optGamma = FALSE
        message('only one rate class, ignored optGamma')
    }
    shape = object$shape
    w = object$w
    g = object$g
    if (type == "DNA" & !is.null(model)) {
        tmp = subsChoice(model)
        optQ = tmp$optQ
        if (!optQ) 
            Q = rep(1, 6)
        optBf = tmp$optBf
        if (!optBf) 
            bf = c(0.25, 0.25, 0.25, 0.25)
        subs = tmp$subs
    }   
    ll0 <- object$logLik
    INV <- object$INV
    ll.0 <- object$ll.0
    rate <- object$rate
    ll = ll0
    ll1 = ll0
    opti = TRUE

    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label))
 
#    on.exit(.C("ll_free"))
#    .C("ll_init", nr, nTips, nc, as.integer(k))
    .INV <- .iind <- NULL
    on.exit({
        pml.free()
#        rm(.INV, .iind)
        })
    pml.init(data, k)    
    
    if (optEdge) {
         res <- optimEdge(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV,
              control = pml.control(epsilon = 1e-07, maxit = 5, trace=trace - 1)) 
         if(trace > 0) 
             cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")  
        if (res[[2]] > ll){  
           ll <- res[[2]]
           tree <- res[[1]]
        }
    }
    if(optRooted){
	    res <- optimRooted(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV, control = pml.control(epsilon = 1e-07, maxit = 10, trace = trace-1))
	    if(trace > 0) 
	        cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
	    if(res[[2]] > ll){  
           ll <- res[[2]]
           tree <- res[[1]]
        }     
	}
    rounds = 1
    while (opti) {
        if (optBf) {
            res = optimBf(tree, data, bf = bf, inv = inv, Q = Q, 
                w = w, g = g, INV = INV, rate = rate, k = k, 
                llMix = llMix)
            bf = res[[1]]
            eig = edQt(Q = Q, bf = bf)
            if (inv > 0) 
                ll.0 <- as.matrix(INV %*% (bf * inv))
            if (wMix > 0) 
                ll.0 <- ll.0 + llMix
            if (trace > 0) 
                cat("optimize base frequencies: ", ll, "-->", 
                  res[[2]], "\n")
            ll = res[[2]]
        }
        if (optQ) {
            if(type=="CODON"){
                 if(is.null(model)) model <- "codon1"
                 model <- match.arg(model, c("codon0", "codon1", "codon2", "codon3"))
                 ab <- c(tstv, dnds)
                 res <- switch(model, 
                     codon1 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab),
                         optK=TRUE, optW = TRUE),  
                     codon2 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab), 
                         optK=FALSE, optW = TRUE),
                     codon3 = optimCodon(tree,data, Q=rep(1,1830), subs=.sub, syn=.syn, 
                         bf = bf, w = w, g = g, inv = inv, INV = INV, ll.0 = ll.0, rate = rate, k = k, ab=log(ab),
                         optK=TRUE, optW = FALSE))
                 tmp <- res[[5]]
                 m = length(tmp)
                 dnds = tmp[m]
                   
                 if(m>1) tstv <- tmp[1]
            }
            else
            res = optimQ(tree, data, Q = Q, subs = subs, bf = bf, w = w, g = g, inv = inv, INV = INV, 
                ll.0 = ll.0, rate = rate, k = k)
            Q = res[[1]]
            eig = edQt(Q = Q, bf = bf)
            if (trace > 0) 
                cat("optimize rate matrix: ", ll, "-->", res[[2]], 
                  "\n")
            ll = res[[2]]
        }
        if(optInv) {
            res = optimInv(tree, data, inv = inv, INV = INV, Q = Q, 
                bf = bf, eig = eig, k = k, shape = shape, rate = rate)
            inv = res[[1]]
            w = rep(1/k, k)
            g = discrete.gamma(shape, k)
            w = (1 - inv) * w
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g = g/(1 - inv)
            g <- g * rate
            ll.0 = as.matrix(INV %*% (bf * inv))
            if (wMix > 0) 
                ll.0 <- ll.0 + llMix
            if (trace > 0) 
                cat("optimize invariant sites: ", ll, "-->", res[[2]], "\n")
            ll = res[[2]]
        }
        if(optGamma) {
            res = optimGamma(tree, data, shape = shape, k = k, 
                inv = inv, INV = INV, Q = Q, bf = bf, eig = eig, 
                ll.0 = ll.0, rate = rate)
            shape = res[[1]]
            w = rep(1/k, k)
            g = discrete.gamma(shape, k)
            if (inv > 0) {
                w = (1 - inv) * w
                g = g/(1 - inv)
            }
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g <- g * rate
            if (trace > 0) 
                cat("optimize shape parameter: ", ll, "-->", 
                  res[[2]], "\n")
            ll = res[[2]]
        }
        if(optRate) {
            res = optimRate(tree, data, rate = rate, inv = inv, 
                INV = INV, Q = Q, bf = bf, eig = eig, k = k, 
                shape = shape, w = w, ll.0 = ll.0)
            if (res[[2]] > ll)rate = res[[1]]
            g = discrete.gamma(shape, k)
            w = rep(1/k, k)
            if (inv > 0) {
                w = (1 - inv) * w
                g = g/(1 - inv)
            }
            if (wMix > 0) 
                w <- (1 - wMix) * w
            g <- g * rate
            if (trace > 0) 
                cat("optimize rate: ", ll, "-->", res[[2]], "\n")
            ll = res[[2]]
        }
        if (optEdge) {  
           res <- optimEdge(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                 control = pml.control(epsilon = 1e-08, maxit = 5, trace=trace - 1)) 
           if (trace > 0) 
              cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
           if (res[[2]] > ll){  
              ll <- res[[2]]
              tree <- res[[1]]
           }
        }
        if(optRooted){
	        res <- optimRooted(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV, control = pml.control(epsilon = 1e-07, maxit = 10, trace = trace-1))
	        if(trace > 0) 
	            cat("optimize edge weights: ", ll, "-->", res[[2]], "\n")
	        if (res[[2]] > ll){  
                ll <- res[[2]]
                tree <- res[[1]]
            }     
	    }
        if(optNni) {
            swap = 0
            iter = 1
            while (iter < 4) {
                if(optEdge){
                    tmp <- pml.nni(tree, data, w, g, eig, bf, ll.0, ll, ...) 
                    swap = swap + tmp$swap
                    res <- optimEdge(tmp$tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, control = pml.control(epsilon = 1e-08, maxit = 3, trace=0)) 
                    ll2 = res[[2]] 
                    tree <- res[[1]]
                }
                else{ 
                    tmp <- rooted.nni(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV, ...) 
                    swap = swap + tmp$swap
                    res <- optimRooted(tmp$tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV, control = pml.control(epsilon = 1e-07, maxit = 5, trace = trace-1))
                    tree <- tmp$tree
                    ll2 = tmp$logLik
                }
                if (trace > 0) 
                  cat("optimize topology: ", ll, "-->", ll2, "\n")
                ll = ll2
                iter = iter + 1
                if (tmp$swap == 0) {
                  iter = 4
                }
            }
            if (trace > 0) 
                cat(swap, "\n")
            if (swap > 0) 
                rounds = 1
            if (swap == 0) 
                optNni = FALSE
        }
        rounds = rounds + 1
        if(rounds > control$maxit) opti <- FALSE
        if (( ll1 - ll ) / ll  < control$eps) #abs(ll1 - ll)
            opti <- FALSE
        ll1 = ll
    }  
    if(type=="CODON"){
        object$dnds = dnds
        object$tstv = tstv
    }
    
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q, 
        levels = attr(data, "levels"), inv = inv, rate = rate, 
        g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, 
        wMix = wMix, site = TRUE)
    
    df <- ifelse(optRooted, tree$Nnode, length(tree$edge.length))
    # length(tree$edge.length)    
    if (type == "CODON") {
        df <- df + (k > 1) + (inv > 0) + 
            length(unique(bf)) - 1 + (dnds != 1) + (tstv != 1) 
    }
    else df = df + (k > 1) + (inv > 0) + 
        length(unique(bf)) - 1 + length(unique(Q)) - 1
    
    object = list(logLik = tmp$loglik, inv = inv, k = k, shape = shape, 
        Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik, weight = attr(data, "weight"), 
        g = g, w = w, eig = eig, data = data, model = model, 
        INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll, 
        call = call, df = df, wMix = wMix, llMix = llMix)
    if (type == "CODON") {
        object$dnds <- dnds
        object$tstv <- tstv
    }
    class(object) = "pml"

    extras = pairlist(bf = bf, Q = Q, inv = inv, shape = shape, rate = rate)[c(optBf, optQ, optInv, optGamma, optRate)]
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    object$call = call   
    object
}


fs <- function (old.el, eig, parent.dat, child.dat, weight, g=g, 
    w=w, bf=bf, ll.0=ll.0, evi, getA=TRUE, getB=TRUE) 
{
    if (old.el < 1e-8) old.el <- 1e-8
    lg = length(parent.dat)
    P <- getP(old.el, eig, g)
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    eve = eig[[2]]
    dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc) 
    X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc) 
    .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), 
            as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)), 
            as.integer(length(weight)), as.double(bf), as.double(weight), 
            as.double(ll.0), as.integer(getA), as.integer(getB))
}


fs3 <- function (old.el, eig, parent.dat, child, weight, g=g, 
    w=w, bf=bf, ll.0=ll.0, contrast, contrast2, evi, ext=TRUE, getA=TRUE, getB=TRUE) # child.dat
{
    if (old.el < 1e-8) old.el <- 1e-8
    lg = length(parent.dat)
    P <- getP(old.el, eig, g)
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    if(ext==FALSE){ 
       child.dat <- child
       eve = eig[[2]]
       dad <- .Call("getDAD", parent.dat, child.dat, P, nr, nc) 
       X <- .Call("getPrep", dad, child.dat, eig[[2]], evi, nr, nc) 
    }
    else {
        nco = as.integer(nrow(contrast))
        dad <- .Call("getDAD2", parent.dat, child, contrast, P, nr, nc, nco)
        child.dat <- vector("list", lg)
        for (i in 1:lg)child.dat[[i]] <- contrast[child, , drop=FALSE]
        X <- .Call("getPrep2", dad, child, contrast2, evi, nr, nc, nco)
    }
    .Call("FS4", eig, as.integer(length(bf)), as.double(old.el), 
            as.double(w), as.double(g), X, child.dat, dad, as.integer(length(w)), 
            as.integer(length(weight)), as.double(bf), as.double(weight), 
            as.double(ll.0), as.integer(getA), as.integer(getB))
}



optimEdge <- function (tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0,
                        control = pml.control(epsilon = 1e-08, maxit = 10, trace=0), ...) 
{
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorder(tree, "postorder") 
    nTips <- length(tree$tip)
    el <- tree$edge.length
    tree$edge.length[el < 1e-08] <- 1e-08
    oldtree = tree
    k = length(w)    
    data = subset(data, tree$tip) 
    loglik = pml.fit2(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
    start.ll <- old.ll <- loglik 
    contrast <- attr(data, "contrast")
    contrast2 <- contrast %*% eig[[2]] 
    evi = (t(eig[[3]]) * bf)
    weight <- attr(data, "weight")
    eps = 1
    iter = 0
    
    treeP = tree
    tree = reorder(tree)
    
    child = tree$edge[, 2]
    parent = tree$edge[, 1]
    m <- max(tree$edge)
    pvec <- integer(m)
    pvec[child] <- parent
    
    EL = numeric(m)
    EL[child] = tree$edge.length
    
    n = length(tree$edge.length)  
    
    nr = as.integer(length(weight))
    nc = as.integer(length(bf))
    nco = as.integer(nrow(contrast))
    eve = eig[[2]]
    lg = k
    rootNode = getRoot(tree)         
    ScaleEPS = 1.0/4294967296.0
    anc = Ancestors(tree, 1:m, "parent")  
    anc0 = as.integer(c(0L, anc))
    
    while (eps > control$eps && iter < control$maxit) {
        blub3 <- .Call("extractScale", as.integer(rootNode), w, g, as.integer(nr), as.integer(nc), as.integer(nTips))
        rowM = apply(blub3, 1, min)       
        blub3 = (blub3-rowM) 
        blub3 = ScaleEPS ^ (blub3) 
        EL <- .Call("optE", as.integer(parent), as.integer(child), 
                    as.integer(anc0), eig, evi, EL, w, g, as.integer(nr), as.integer(nc), 
                    as.integer(nTips), as.double(contrast), 
                    as.double(contrast2), nco, blub3, data, as.double(weight), as.double(ll.0))       
        iter = iter + 1
#        tree$edge.length = EL[tree$edge[,2]]
        treeP$edge.length = EL[treeP$edge[,2]]
        newll <- pml.fit2(treeP, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k)
        
        eps = ( old.ll - newll ) / newll
        if( eps <0 ) return(list(oldtree, old.ll))
        oldtree = treeP
        if(control$trace>1) cat(old.ll, " -> ", newll, "\n") 
        old.ll = newll
        #        loli = parent[1] 
    }
    if(control$trace>0) cat(start.ll, " -> ", newll, "\n")
    list(tree=treeP, logLik=newll, c(eps, iter))
}


# bf raus C naeher
# data=data, k=k, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV)
pml.move <- function(EDGE, el, data, g, w, eig, k, nTips, bf){
    node <- EDGE[, 1]
    edge <- EDGE[, 2]
    root <- as.integer(node[length(node)])     
#    el <- as.double(tree$edge.length)
    nr = as.integer(attr(data, "nr"))
    nc = as.integer(attr(data, "nc"))    
    node = as.integer(node - nTips - 1L)  
    edge = as.integer(edge - 1L)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])    
    tmp <- .Call("PML3", dlist=data, as.double(el), as.double(w), as.double(g), nr, nc, k, eig, as.double(bf), node, edge, nTips, root, nco, contrast, N=as.integer(length(edge))) 
    return(NULL)
}



#
# pmlPart + pmlCluster
#
optimPartQ <- function (object, Q = c(1, 1, 1, 1, 1, 1), ...) 
{
    l = length(Q)
    Q = Q[-l]
    Q = sqrt(Q)
    fn = function(Q, object, ...) {
        result <- 0
        Q = c(Q^2, 1)
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], Q = Q, ...)$logLik
        result
    }
    res = optim(par = Q, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = 0, upper = Inf, control = list(fnscale = -1, 
            maxit = 25), object = object, ...)
    res[[1]] = c(res[[1]]^2, 1)
    res
}


optimPartQGeneral <- function (object, Q = c(1, 1, 1, 1, 1, 1), subs=rep(1,length(Q)), ...) 
{
    m = length(Q)
    n = max(subs)
    ab = numeric(n)
    for(i in 1:n) ab[i]=log(Q[which(subs==i)[1]])
    fn = function(ab, object, m, n, subs, ...) {
        Q = numeric(m)
        for(i in 1:n)Q[subs==i] = ab[i]
        Q = exp(Q)
        result = 0
        for (i in 1:length(object)) result <- result + update(object[[i]], Q = Q, ...)$logLik
        result
    }
    res = optim(par = ab, fn = fn, gr = NULL, method = "L-BFGS-B", 
        lower = -Inf, upper = Inf, control = list(fnscale = -1, 
            maxit = 25), object = object, m=m, n=n, subs=subs, ...)
    Q = rep(1, m)
    for(i in 1:n) Q[subs==i] = exp(res[[1]][i])
    res[[1]] = Q
    res
}


optimPartBf <- function (object, bf = c(0.25, 0.25, 0.25, 0.25), ...) 
{
    l = length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, object, ...) {
        result <- 0
        bf = exp(c(lbf, 0))
        bf = bf/sum(bf)
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], 
            bf = bf, ...)$logLik
        result
    }
    res = optim(par = lbf, fn = fn, gr = NULL, method = "Nelder-Mead", 
        control = list(fnscale = -1, maxit = 500), object, ...)
    print(res[[2]])
    bf = exp(c(res[[1]], 0))
    bf = bf/sum(bf)
}


optimPartInv <- function (object, inv = 0.01, ...) 
{
    fn = function(inv, object, ...) {
        result <- 0
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], inv = inv, 
            ...)$logLik
        result
    }
    res = optimize(f = fn, interval = c(0, 1), lower = 0, upper = 1, 
        maximum = TRUE, tol = 1e-04, object, ...)
#    print(res[[2]])
    res[[1]]
}


optimPartGamma <- function (object, shape = 1, ...) 
{
    fn = function(shape, object, ...) {
        result <- 0
        n <- length(object)
        for (i in 1:n) result <- result + update(object[[i]], shape = shape, 
            ...)$logLik
        result
    }    
    res = optimize(f = fn, interval = c(0, 100), lower = 0, upper = 100, 
        maximum = TRUE, tol = 0.01, object, ...)
    res
}


dltmp <- function (fit, i=1, transform=transform) # i = weights
{
    tree = fit$tree 
    data = getCols(fit$data, tree$tip) 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    q = length(tree$tip.label)
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    m = length(edge) + 1  # max(edge)
    dat = vector(mode = "list", length = m)
    eig = fit$eig
    w = fit$w[i]
    g = fit$g[i]
    bf = fit$bf
    el <- tree$edge.length
    P <- getP(el, eig, g)
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
    node = as.integer(node - min(node))
    edge = as.integer(edge - 1)
    nTips = as.integer(length(tree$tip))
    mNodes = as.integer(max(node) + 1)
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])
    dat[(q + 1):m] <- .Call("LogLik2", data, P, nr, nc, node, edge, nTips, mNodes, contrast, nco)
    result = dat[[q+1]] %*% (bf * w)

    parent <- tree$edge[, 1]
    child <- tree$edge[, 2]
    nTips = min(parent) - 1
    datp = vector("list", m)
    el = tree$edge.length 
    if (transform) dP = getdP(tree$edge.length, eig, g)
    else dP = getdP2(tree$edge.length, eig, g)
   
    datp[(nTips + 1)] = dat[(nTips + 1)]
    l = length(child)
    dl = matrix(0, nr, l)
    for (j in (m - 1):1) {
        # tips have factor format, internal edges are matrices
        if (child[j] > nTips){
             tmp2 = (datp[[parent[j]]]/(dat[[child[j]]] %*% P[[j]]))
             dl[, j] = (tmp2 * (dat[[child[j]]] %*% dP[[j]])) %*% (w * bf)
             datp[[child[j]]] = (tmp2 %*% P[[j]]) * dat[[child[j]]]  
             }
        else{
             tmp2 = (datp[[parent[j]]]/((contrast %*% P[[j]])[data[[child[j]]],] ))
             dl[, j] = (tmp2 * ((contrast %*% dP[[j]])[data[[child[j]]],]) ) %*% (w * bf)    
             }
    }
    dl
}


dl <- function(x, transform = TRUE){
  w = x$w 
  l=length(x$w)
  dl = dltmp(x, 1, transform)
  i=2
  while(i < (l+1)){
    dl = dl + dltmp(x, i, transform)
    i = i + 1
  } 
  dl
}


# add control and change edge
optimPartEdge <- function (object, ...) 
{
    tree <- object[[1]]$tree
    theta <- object[[1]]$tree$edge.length
    n <- length(object)
    l <- length(theta)
    nrv <- numeric(n)
    for (i in 1:n) nrv[i] = attr(object[[i]]$data, "nr")
    cnr <- cumsum(c(0, nrv))
    weight = numeric(sum(nrv))
    dl <- matrix(NA, sum(nrv), l)
    for (i in 1:n) weight[(cnr[i] + 1):cnr[i + 1]] = attr(object[[i]]$data, 
        "weight")
    ll0 = 0
    for (i in 1:n) ll0 = ll0 + object[[i]]$logLik
    eps = 1
    scalep =1
    k = 1
    while (eps > 0.001 & k<50) {
        if(scalep==1){
            for (i in 1:n) {
                lv = drop(exp(object[[i]]$site))
                dl[(cnr[i] + 1):cnr[i + 1], ] = dl(object[[i]], TRUE)/lv
            }
            sc = colSums(weight * dl)
            F = crossprod(dl * weight, dl) + diag(l)*1e-10
            # add small ridge penalty for numerical stability 
        }
        thetaNew = log(theta) + scalep * solve(F, sc)
        tree$edge.length = as.numeric(exp(thetaNew))
        for (i in 1:n) object[[i]] <- update(object[[i]], tree = tree)
        ll1 = 0
        for (i in 1:n) ll1 = ll1 + object[[i]]$logLik
        eps <- ll1 - ll0
        if (eps < 0 || is.nan(eps)) {
            scalep = scalep/2
            eps = 1
            thetaNew = log(theta)
            ll1 = ll0
        }
        else scalep = 1
        theta = exp(thetaNew)
        ll0 <- ll1
        k=k+1
    }
    object
}


makePart <- function(fit, rooted, weight=~index+genes){
    if(class(fit)=="phyDat"){
        x <- fit
        dm <- dist.ml(x)
        if(!rooted) tree <- NJ(dm)
        else tree <- upgma(dm)
        fit <- pml(tree, x, k=4)
    }     
    dat <- fit$data 
    if(class(weight)[1]=="formula")     
        weight <- xtabs(weight, data=attr(dat, "index"))
    fits <- NULL 
    for(i in 1:dim(weight)[2]){ 
        ind <- which(weight[,i] > 0)
        dat2 <- getRows(dat, ind)
        attr(dat2, "weight") <- weight[ind,i]
        fits[[i]] <- update(fit, data = dat2)
    }
    names(fits) = colnames(fits)
    fits    
}


multiphyDat2pmlPart <- function(x, rooted=FALSE, ...){
    fun <-  function(x, ...){
        dm <- dist.ml(x)
        if(!rooted) tree <- NJ(dm)
        else tree <- upgma(dm)
        fit <- pml(tree, x, ...)
    }
    fits <- lapply(x@dna, fun, ...)
    fits
}


pmlPart2multiPhylo <- function(x){
    res <- lapply(x$fits, FUN=function(x)x$tree)
    class(res) <- "multiPhylo"
    res
}


plot.pmlPart<- function(x, ...){
    plot(pmlPart2multiPhylo(x), ...)
}


pmlPart <- function (formula, object, control=pml.control(epsilon=1e-8, maxit=10, trace=1), model=NULL, rooted=FALSE, ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllNNI <- !is.na(optAll[1])
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    PartNni <- !is.na(optPart[1])
    PartBf <- !is.na(optPart[2])
    PartQ <- !is.na(optPart[3])
    PartInv <- !is.na(optPart[4])
    PartGamma <- !is.na(optPart[5])
    PartEdge <- !is.na(optPart[6])
    PartRate <- !is.na(optPart[7])
 
    if(class(object)=="multiphyDat"){
        if(AllNNI || AllEdge) object <- do.call(cbind.phyDat, object@dna)
        else fits <- multiphyDat2pmlPart(object, rooted=rooted, ...)
    } 
    if(class(object)=="pml") fits <- makePart(object, rooted=rooted, ...) 
    if(class(object)=="phyDat") fits <- makePart(object, rooted=rooted, ...)
    if(class(object)=="pmlPart") fits <- object$fits
    if(class(object)=="list") fits <- object


    trace = control$trace
    epsilon = control$epsilon
    maxit = control$maxit

    p <- length(fits)
 #   if(length(model)<p) model = rep(model, length = p)

    m = 1
    logLik = 0
    for (i in 1:p) logLik = logLik + fits[[i]]$log
    eps = 10
    while (eps > epsilon & m < maxit) {
        loli = 0
        if(any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
            for (i in 1:p) {
                fits[[i]] = optim.pml(fits[[i]], optNni=PartNni, optBf=PartBf, 
                    optQ=PartQ, optInv=PartInv, optGamma=PartGamma,  optEdge=PartEdge, 
                    optRate=PartRate, optRooted=rooted,  
                    control = pml.control(maxit = 3, epsilon = 1e-8, trace-1), model=model[i])
            }
        } 
        if (AllQ) {
            Q = fits[[1]]$Q
            subs = c(1:(length(Q)-1), 0)
            newQ <- optimPartQGeneral(fits, Q=Q, subs=subs)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ[[1]])
        }
        if (AllBf) {
             bf = fits[[1]]$bf
            newBf <- optimPartBf(fits, bf=bf)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
        }
        if (AllInv) {
            inv = fits[[1]]$inv
            newInv <- optimPartInv(fits, inv=inv)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv)
        }
        if (AllGamma) {
            shape = fits[[1]]$shape
            newGamma <- optimPartGamma(fits, shape=shape)[[1]]
            for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
        }
        if (AllNNI){
            fits <- optimPartNNI(fits,AllEdge)
            if(trace>0) cat(attr(fits,"swap"), " NNI operations performed")
        }
        if (AllEdge) 
            fits <- optimPartEdge(fits)
        if (PartRate){
            tree = fits[[1]]$tree
            rate=numeric(p)
            wp =numeric(p) 
            for(i in 1:p){
                wp[i]=sum(fits[[i]]$weight)
                rate[i] <- fits[[i]]$rate
                }          
            ratemult = sum(wp) / sum(wp*rate)
            tree$edge.length = tree$edge.length/ratemult  
            for(i in 1:p)fits[[i]] = update(fits[[i]], tree=tree, rate=rate[i]*ratemult)   
        }
        loli <- 0
        for (i in 1:p) loli <- loli + fits[[i]]$log
        eps = (logLik - loli)/loli
        if(trace>0) cat("loglik:", logLik, "-->", loli, "\n")
        logLik <- loli
        m = m + 1
    }
    
    df <- matrix(1, 6 ,2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1,1] <- length(fits[[1]]$tree$edge.length)
    df[2,1] <- fits[[1]]$k > 1
    df[3,1] <- fits[[1]]$inv > 0
    df[4,1] <- length(unique(fits[[1]]$bf)) - 1
    df[5,1] <- length(unique(fits[[1]]$Q)) - 1
    df[6,1] <- 0 # rates 
    if(PartEdge) df[1,2] = p
    if(PartGamma) df[2,2] = p
    if(PartInv) df[3,2] = p
    if(PartBf) df[4,2] = p
    if(PartQ) df[5,2] = p
    if(PartRate) df[6,1] = p-1     
    attr(logLik, "df") = sum(df[,1]*df[,2])
    object <- list(logLik = logLik, fits = fits, call = call, df=df)
    class(object) <- "pmlPart" 
    object
}



bip <- function (obj) 
{
    if (is.null(attr(obj, "order")) || attr(obj, "order") == 
            "cladewise") 
        obj <- reorder(obj, "postorder")
    maxP = max(obj$edge)
    nTips = length(obj$tip)
    res <- .Call("C_bip", as.integer(obj$edge[, 1]), as.integer(obj$edge[, 2]), as.integer(nTips), as.integer(maxP))
    res
}


bipart <- function(obj){
    if (is.null(attr(obj, "order")) || attr(obj, "order") == "cladewise") 
        obj <- reorder(obj, "postorder")
    maxP  = max(obj$edge)
    nTips = length(obj$tip)
    res <- .Call("C_bipart", as.integer(obj$edge[,1]) , as.integer(obj$edge[,2]), as.integer(nTips), as.integer(maxP))  #, as.integer(obj$Nnode))
#    attr(res, "nodes") = unique(obj$edge[,1])
    res    
}


bipartition <- function (tree) 
{
    if(is.rooted(tree))tree <- unroot(tree)
    if(is.null(attr(tree,"order")) || attr(tree, "order")=="cladewise") tree <- reorder(tree, "postorder")
    bp <- bipart(tree)
    nTips = length(tree$tip)
    l = length(bp)
    m = length(bp[[l]])
    k = length(tree$edge[, 1])
    result = matrix(0L, l, m)
    res = matrix(0L, k, m)
    for (i in 1:l) result[i, bp[[i]]] = 1L
    result = result[-l, ,drop=FALSE]
    for (i in 1:nTips) res[(tree$edge[, 2] == i), i] = 1L     
#    res[tree$edge[, 2] > nTips, ] = result
    res[ match(unique(tree$edge[,1]),tree$edge[,2])[-l], ] = result
    colnames(res) = tree$tip.label
    rownames(res) = tree$edge[,2]
    res[res[, 1] == 1, ] = 1L - res[res[, 1] == 1, ]
    res
}



pmlCluster.fit <- function (formula, fit, weight, p = 4, part = NULL, control=pml.control(epsilon=1e-8, maxit=10, trace=1), ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllNNI <- !is.na(optAll[1])
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    PartNni <- !is.na(optPart[1])
    PartBf <- !is.na(optPart[2])
    PartQ <- !is.na(optPart[3])
    PartInv <- !is.na(optPart[4])
    PartGamma <- !is.na(optPart[5])
    PartEdge <- !is.na(optPart[6])
    PartRate <- !is.na(optPart[7])
    nrw <- dim(weight)[1]
    ncw <- dim(weight)[2]
    if (is.null(part)){ 
        part = rep(1:p, length=ncw)
        part = sample(part)
        }
    Part = part
    Gtrees = vector("list", p)
    dat <- fit$data
    attr(fit$orig.data, "index") <- attr(dat, "index") <- NULL
    for (i in 1:p) Gtrees[[i]] = fit$tree
    fits = vector("list", p)
    for (i in 1:p) fits[[i]] = fit
    trace = control$trace
    eps = 0
    m = 1
    logLik = fit$log
    trees = list()
    weights = matrix(0, nrw, p)
    lls = matrix(0, nrw, p)
    loli = fit$log
    oldpart = part
    eps2 = 1
    iter = 0
    swap = 1
    while (eps < ncw || abs(eps2) > control$eps) {
        df2 = 0
        
        if(any(c(PartNni, PartBf, PartInv, PartQ, PartGamma, PartEdge, PartRate))){
            for (i in 1:p) {
                weights[, i] = rowSums(weight[, which(part == i), 
                    drop = FALSE])
                ind <- which(weights[, i] > 0)
                dat2 <- getRows(dat, ind)
                attr(dat2, "weight") <- weights[ind, i]
                fits[[i]] <- update(fits[[i]], data = dat2)
                fits[[i]] = optim.pml(fits[[i]], PartNni, PartBf, 
                    PartQ, PartInv, PartGamma, PartEdge, PartRate, 
                    control = pml.control(epsilon = 1e-8, maxit = 3, trace-1))
                lls[, i] = update(fits[[i]], data = dat)$site
                Gtrees[[i]] = fits[[i]]$tree
            }
        }
        if (AllQ) {
            Q = fits[[1]]$Q
            subs = c(1:(length(Q)-1), 0)
            newQ <- optimPartQGeneral(fits, Q=Q, subs=subs)[[1]]
            for (i in 1:p) fits[[i]] <- update(fits[[i]], Q = newQ)
            df2 = df2 + length(unique(newQ)) - 1
        }
        if (AllBf) {
	        bf = fits[[1]]$bf
            newBf <- optimPartBf(fits, bf=bf)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], bf = newBf)
            df2 = df2 + length(unique(newBf)) - 1
        }
        if (AllInv) {
            inv = fits[[1]]$inv
            newInv <- optimPartInv(fits, inv=inv)
            for (i in 1:p) fits[[i]] <- update(fits[[i]], inv = newInv) #there was an Error
            df2 = df2 + 1
        }
        if (AllGamma) {
            shape = fits[[1]]$shape
            newGamma <- optimPartGamma(fits, shape=shape)[[1]]        
            for (i in 1:p) fits[[i]] <- update(fits[[i]], shape = newGamma)
            df2 = df2 + 1
        }
        if (AllNNI) {
            fits <- optimPartNNI(fits, AllEdge)
            if(trace>0)cat(attr(fits, "swap"), " NNI operations performed")
            swap <- attr(fits, "swap")
        }
        if (AllEdge) {
            fits <- optimPartEdge(fits)
            df2 = df2 + length(fits[[1]]$tree$edge.length)
        }
        if (PartRate) {
            tree = fits[[1]]$tree
            rate = numeric(p)
            wp = numeric(p)
            for (i in 1:p) {
                wp[i] = sum(fits[[i]]$weight)
                rate[i] <- fits[[i]]$rate
            }
            ratemult = sum(wp)/sum(wp * rate)
            tree$edge.length = tree$edge.length/ratemult
            for (i in 1:p) fits[[i]] = update(fits[[i]], tree = tree, 
                rate = rate[i] * ratemult)
        }
        for (i in 1:p) lls[, i] = update(fits[[i]], data = dat)$site
        trees[[m]] = Gtrees
        LL = t(weight) %*% lls       
# choose partitions which change        
        tmp =(LL[cbind(1:ncw,part)] - apply(LL, 1, max))/colSums(weight)
        fixi = numeric(p)
        for(i in 1:p){
            tmpi = which(part == i)
            fixi[i] = tmpi[which.max(tmp[tmpi])]     
            }
        oldpart = part
# restrict the number of elements changing groups 
# If more than 25% would change, only the 25% with the highest increase per site change       
        if( sum(tmp==0)/length(tmp) < .75){
           medtmp = quantile(tmp, .25)
           medind = which(tmp<=medtmp)
           part[medind] = apply(LL[medind,], 1, which.max)
           }
        else part = apply(LL, 1, which.max)
# force groups to have at least one member
        part[fixi] = 1:p
        Part = cbind(Part, part)
        eps = sum(diag(table(part, oldpart)))
        eps2 = loli
        loli = sum(apply(LL, 1, max))
        eps2 = (eps2 - loli)/loli
        logLik = c(logLik, loli)
        if(trace>0) print(loli)
        Part = cbind(Part, part)
        df2 = df2 + df2
        if (eps == ncw & swap == 0) 
            AllNNI = FALSE
        m = m + 1
        if (eps == ncw) 
            iter = iter + 1
        if (iter == 3) 
            break
    }
    df <- matrix(1, 6, 2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1, 1] <- length(fits[[1]]$tree$edge.length)
    df[2, 1] <- fits[[1]]$k - 1
    df[3, 1] <- fits[[1]]$inv > 0
    df[4, 1] <- length(unique(fits[[1]]$bf)) - 1
    df[5, 1] <- length(unique(fits[[1]]$Q)) - 1
    df[6, 1] <- 0
    if (PartEdge) 
        df[1, 2] = p
    if (PartGamma) 
        df[2, 2] = p
    if (PartInv) 
        df[3, 2] = p
    if (PartBf) 
        df[4, 2] = p
    if (PartQ) 
        df[5, 2] = p
    if (PartRate) 
        df[6, 1] = p - 1
    attr(logLik, "df") = sum(df[, 1] * df[, 2])
    res = list(logLik = logLik, Partition = Part, trees = trees) # intermediate results
    result <- list(logLik = loli, fits = fits, Partition = part, df = df, res = res, call = call)
    class(result) <- c("pmlPart")
    result
}


pmlCluster <- function (formula, fit, weight, p = 1:5, part = NULL, nrep = 10, control = pml.control(epsilon = 1e-08,
   maxit = 10, trace = 1), ...)
{
   call <- match.call()
   form <- phangornParseFormula(formula)
   if(any(p==1)){
       opt2 <- c("nni", "bf", "Q", "inv", "shape", "edge")
       tmp1 <- opt2 %in% form$left
       tmp1 <- tmp1 | (opt2 %in% form$right)
       fit <- optim.pml(fit, tmp1[1], tmp1[2], tmp1[3], tmp1[4],
       tmp1[5], tmp1[6])
   }

   p=p[p!=1]
   if(length(p)==0)return(fit)
   n = sum(weight)
   k=2

   BIC = matrix(0, length(p)+1, nrep)
   BIC[1,] = AIC(fit, k = log(n))
   LL = matrix(NA, length(p)+1, nrep)
   LL[1,] = logLik(fit)

   P = array(dim=c(length(p)+1, nrep, dim(weight)[2]))
   tmpBIC = Inf
   choice = c(1,1) 
   for(j in p){
       tmp=NULL
       for(i in 1:nrep){
           tmp = pmlCluster.fit(formula, fit, weight, p=j, part=part, control=control,...)
           P[k,i,] = tmp$Partition
           BIC[k,i] = AIC(tmp, k = log(n))
           LL[k,i] = logLik(tmp)
           if(BIC[k,i]<tmpBIC){
                tmpBIC = BIC[k,i]
                result = tmp
                choice = c(k,i) 
           }
       }
       k=k+1
   }      

   p = c(1,p)
   result$choice = choice 
   result$BIC = BIC
   result$AllPartitions = P
   result$AllLL = LL
   result$p = p 
   class(result) = c("pmlCluster", "pmlPart")
   result
}


plot.pmlCluster <- function(x, which = c(1L:3L), caption = list("BIC", "log-likelihood", "Partitions"), ...){
   show <- rep(FALSE, 3)
   show[which] <- TRUE
   choice = x$choice
   if(show[1]){
       X <- x$AllPartitions[choice[1],,]
       d <- dim(X)
       ind = order(X[choice[2],])
       im  = matrix(0, d[2], d[2])
       for(j in 1:d[1]){for(i in 1:d[2]) im[i,] = im[i,] + (X[j,] == X[j,i]) }
       image(im[ind, ind], ...)
   }

   if(show[1])matplot(x$p, x$BIC, ylab="BIC", xlab="number of clusters")
   if(show[1])matplot(x$p, x$AllLL, ylab="log-likelihood", xlab="number of clusters")
}


readAArate <- function(file){
    tmp <- read.table(system.file(file.path("extdata", file)), col.names = 1:20, fill=TRUE)
    Q <- tmp[1:19,1:19]
    names <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l", "k", "m", "f", "p", "s", "t", "w",  "y", "v")
    Q <- as.numeric(Q[lower.tri(Q,TRUE)])
    bf <- as.numeric(as.character(unlist(tmp[20,])))
    names(bf) <- names
    list(Q=Q, bf=bf)
}


#.LG <- readAArate("lg.dat")
#.WAG <- readAArate("wag.dat")
#.Dayhoff <- readAArate("dayhoff-dcmut.dat")
#.JTT <- readAArate("jtt-dcmut.dat")
#.cpREV <- readAArate("cpREV.dat")
#.mtmam <- readAArate("mtmam.dat")
#.mtArt <- readAArate("mtArt.dat")
# save(.LG,.WAG,.Dayhoff,.JTT,.cpREV,.mtmam,.mtArt, file = "sysdata2.rda")


getModelAA <- function(model, bf=TRUE, Q=TRUE){
    model <- match.arg(eval(model), .aamodels)
    tmp = get(paste(".", model, sep=""), environment(pml))
    if(Q) assign("Q", tmp$Q, envir=parent.frame())
    if(bf) assign("bf", tmp$bf, envir=parent.frame())
}


print.pml = function(x,...){
    cat("\n loglikelihood:", x$logLik, "\n")
    w <- x$weight
    w <- w[w>0]    
    type <- attr(x$data, "type")
    levels <- attr(x$data, "levels")
    nc <- attr(x$data, "nc")
    ll0 = sum(w*log(w/sum(w)))
    cat("\nunconstrained loglikelihood:", ll0, "\n")
    if(x$inv > 0)cat("Proportion of invariant sites:",x$inv,"\n")
    if(x$k >1){
        cat("Discrete gamma model\n")
        cat("Number of rate categories:",x$k,"\n")        
        cat("Shape parameter:",x$shape,"\n")
        }
    if(type=="AA") cat("Rate matrix:",x$model, "\n")    
    if(type=="DNA"){
        cat("\nRate matrix:\n")    
        QM = matrix(0, nc, nc, dimnames = list(levels,levels))    
        QM[lower.tri(QM)] = x$Q    
        QM = QM+t(QM)
        print(QM)
        cat("\nBase frequencies:  \n")
        bf = x$bf
        names(bf) = levels 
        cat(bf, "\n")
    }
    if(type=="CODON") {
         cat("dn/ds:",x$dnds, "\n")
         cat("ts/tv:",x$tstv, "\n") 
    }
    if(type=="USER" & length(x$bf)<11){         
        cat("\nRate matrix:\n")    
        QM = matrix(0, nc, nc, dimnames = list(levels,levels))    
        QM[lower.tri(QM)] = x$Q    
        QM = QM+t(QM)
        print(QM)
        cat("\nBase frequencies:  \n")
        bf = x$bf
        names(bf) = levels 
        cat(bf, "\n")
    }        
}


optEdgeMulti <- function (object, control = pml.control(epsilon = 1e-8, maxit = 10, trace=1), ...) 
{
    tree <- object$tree
    theta <- object$tree$edge.length
    weight <- attr(object$data, "weight")
    ll0 = object$logLik
    eps = 1
    iter = 0
    iter2 = 0
    scale = 1
    # l = length(theta)
    while (abs(eps) > control$eps && iter < control$maxit) {
        dl = score(object)
        thetaNew = log(theta) + scale * solve(dl[[2]], dl[[1]]) #+ diag(l)*1e-10
        newtheta = exp(thetaNew)
        tree$edge.length = as.numeric(newtheta)
        object <- update(object, tree = tree)
        ll1 = object$logLik 
        eps <- ( ll0 - ll1 ) / ll1 
        if(eps < 0){
             newtheta = theta
             scale = scale / 2
             tree$edge.length = as.numeric(theta)  
             ll1 = ll0  
             iter2 <- iter2+1             
        }
        else{
            scale=1
            iter2 = 0
        }  
        theta = newtheta 
        if(iter2==0 && control$trace>0) cat("loglik: ",ll1,"\n")
        ll0 <- ll1
        if(iter2==10)iter2=0  
        if(iter2==0)iter <- iter+1
    }
    object <- update(object, tree = tree) 
    object
}


# add data for internal use parent.frame(n) for higher nestings 
update.pmlNew <- function (object, ..., evaluate = TRUE){
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) 
        eval(call, object, parent.frame())
    else call
}


update.pml <- function (object, ...) 
{
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("tree", "data", "bf", "Q", "inv", "k", "shape", 
        "rate", "model", "wMix", "llMix", "...") 
    names(extras) <- pmla[pmatch(names(extras), pmla[-length(pmla)])]
    call = object$call
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }    
    existing <- match(pmla, names(extras))
    updateEig <- FALSE
    updateRates <- FALSE
    if (is.na(existing[1])) tree <- object$tree
    else tree <- eval(extras[[existing[1]]], parent.frame())
    if(is.null(attr(tree,"order")) || attr(tree,"order")=="cladewise")tree <- reorder(tree, "postorder")
    if (is.na(existing[2])){
        data <- object$data
        INV <- object$INV
        }
    else{ 
        data <- eval(extras[[existing[2]]], parent.frame())
        ll.0 <- numeric(attr(data,"nr"))
        INV <- Matrix(lli(data, tree), sparse=TRUE)
    }
    nr <- as.integer(attr(data, "nr"))
    nc <- as.integer(attr(data, "nc"))
      
    if (is.na(existing[3])) bf <- object$bf
    else {
        bf <- eval(extras[[existing[3]]], parent.frame())
        updateEig <- TRUE
    }
    if (is.na(existing[4])) Q <- object$Q
    else {
         Q <- eval(extras[[existing[4]]], parent.frame())
         updateEig <- TRUE  
    }    
#    model <- object$model
    type <- attr(object$data, "type")
    model<-NULL
    if (type == "AA") {
        if(!is.na(existing[9]) ){  
#        model <- match.arg(eval(extras[[existing[9]]], parent.frame()), c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        model <- match.arg(eval(extras[[existing[9]]], parent.frame()), .aamodels)
        getModelAA(model, bf = is.na(existing[3]), Q = is.na(existing[4]))
        updateEig <- TRUE
        } 
#        else model <- object$model
    }
   
    if(is.na(existing[5])) inv <- object$inv
    else{
        inv <- eval(extras[[existing[5]]], parent.frame())
        updateRates <- TRUE
    }
    if(is.na(existing[6])) k <- object$k
    else{
        k <- eval(extras[[existing[6]]], parent.frame())
        updateRates <- TRUE
    }
    if(is.na(existing[7])) shape <- object$shape
    else{
        shape <- eval(extras[[existing[7]]], parent.frame())
        updateRates <- TRUE
    }
    rate <- ifelse(is.na(existing[8]), object$rate, eval(extras[[existing[8]]], parent.frame()))
    wMix <- ifelse(is.na(existing[10]), object$wMix, eval(extras[[existing[10]]], parent.frame()))
    if(is.na(existing[11])) llMix <- object$llMix
    else llMix <- eval(extras[[existing[11]]], parent.frame())
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    if(updateEig)eig <- edQt(bf = bf, Q = Q)
    else eig <- object$eig
    g <- discrete.gamma(shape, k)
    g <- rate * g 
    if (inv > 0) g <- g/(1 - inv)
    ll.0 <- as.matrix(INV %*% (bf * inv))
    if(wMix>0) ll.0 <- ll.0 + llMix
    w = rep(1/k, k)
    if (inv > 0) 
        w <- (1 - inv) * w
    if (wMix > 0) 
        w <- wMix * w                  
    m <- 1
    
    resll <- matrix(0, nr, k)
    nTips = as.integer(length(tree$tip.label))  

    data <- subset(data, tree$tip.label)     

    on.exit(.C("ll_free"))
    .C("ll_init", nr, nTips, nc, as.integer(k))
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q, levels = attr(data, "levels"),
        inv = inv, rate = rate, g = g, w = w, eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix,
        wMix = wMix, site = TRUE)
    
    df <- ifelse(is.ultrametric(tree), tree$Nnode, length(tree$edge.length))
    
    if (type == "CODON") {
        df <- df + (k > 1) + (inv > 0) + length(unique(bf)) - 1
    }
    else df = df + (k > 1) + (inv > 0) + 
        length(unique(bf)) - 1 + length(unique(Q)) - 1
    result = list(logLik = tmp$loglik, inv = inv, k = k, shape = shape, Q = Q, bf = bf, 
        rate = rate, siteLik = tmp$siteLik, weight = weight, g = g, w = w, eig = eig, 
        data = data, model = model, INV = INV, ll.0 = ll.0, tree = tree, lv = tmp$resll,
        call = call, df = df, wMix = wMix, llMix = llMix)
    if (type == "CODON") {
        result$dnds <- 1
        result$tstv <- 1
    }
    class(result) = "pml"
    result 
}


optimMixQ <- function(object, Q=c(1, 1, 1, 1, 1, 1), omega,...){
    l = length(Q)
    Q = Q[-l]
    Q = sqrt(Q)
    fn = function(Q, object, omega,...) {
        Q = c(Q^2, 1)
        weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + as.numeric(update(object[[i]], Q=Q, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optim(par=Q, fn=fn, gr=NULL, method="L-BFGS-B", lower=0, 
            upper=Inf, control=list(fnscale = -1, maxit=25), 
            object=object, omega=omega,...)
    res[[1]] = c(res[[1]]^2, 1)
    res
}


optimMixBf <- function(object, bf=c(.25,.25,.25,.25), omega,...){
    l = length(bf)
    nenner = 1/bf[l]
    lbf = log(bf * nenner)
    lbf = lbf[-l]
    fn = function(lbf, object, omega,...) {
    bf = exp(c(lbf,0))
    bf = bf/sum(bf)
    weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + as.numeric(update(object[[i]], bf=bf, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optim(par=lbf, fn=fn, gr=NULL, method="Nelder-Mead", 
        control=list(fnscale=-1, maxit=500), object, omega=omega,...)
#    print(res[[2]])
    bf = exp(c(res[[1]],0))
    bf = bf/sum(bf)
}


optimMixInv <- function(object, inv=0.01, omega,...){
    fn = function(inv, object, omega,...) {
        n <- length(omega)
        weight <- object[[1]]$weight
        p <- length(weight)
        result <- numeric(p)
         for(i in 1:n)result <- result + as.numeric(update(object, inv=inv, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res = optimize(f=fn, interval = c(0,1), lower = 0, upper = 1, maximum = TRUE,
        tol = .0001, object, omega=omega,...)
#    print(res[[2]]) 
    res[[1]]
}


optimMixRate <- function (fits, ll, weight, omega, rate=rep(1,length(fits))) 
{
    r <- length(fits)
    rate0 <- rate[-r]   

    fn<-function(rate, fits, ll, weight, omega){
        r <-  length(fits)
        rate <- c(rate, (1- sum(rate *omega[-r]))/omega[r])
        for (i in 1:r) fits[[i]]<-update(fits[[i]], rate = rate[i])
        for (i in 1:r) ll[, i] <- fits[[i]]$lv
        sum(weight*log(ll%*%omega)) 
    }
    ui=diag(r-1)
    ui <- rbind(-omega[-r], ui)
    ci <- c(-1, rep(0, r-1))
    res <- constrOptim(rate0, fn, grad=NULL, ui=ui, ci=ci, mu = 1e-04, control = list(fnscale=-1),
        method = "Nelder-Mead", outer.iterations = 50, outer.eps = 1e-05, fits=fits, ll=ll, weight=weight, omega=omega)
    rate <- res[[1]]
    res[[1]] <- c(rate, (1- sum(rate *omega[-r]))/omega[r])
    res
}


optW <- function (ll, weight, omega,...) 
{
    k = length(omega)
    nenner = 1/omega[1]
    eta = log(omega * nenner)
    eta = eta[-1]
    fn = function(eta, ll, weight) {
        eta = c(0,eta)
        p = exp(eta)/sum(exp(eta))
        res = sum(weight * log(ll %*% p)) 
        res
    }
    if(k==2)res = optimize(f =fn , interval =c(-3,3) , lower = -3, upper = 3, maximum = TRUE, tol = .Machine$double.eps^0.25, ll = ll, weight = weight) 
    else res = optim(eta, fn = fn, method = "L-BFGS-B", lower=-5, upper=5,control = list(fnscale = -1, 
        maxit=25), gr = NULL, ll = ll, weight = weight)

    p = exp(c(0,res[[1]]))
    p = p/sum(p)
    result = list(par = p, value = res[[2]])
    result
}


optimMixEdge <- function(object, omega, trace=1,...){
    tree <- object[[1]]$tree
    theta <- object[[1]]$tree$edge.length
    weight = as.numeric(attr(object[[1]]$data,"weight"))
    n <- length(omega)
    p <- length(weight)
    q <- length(theta)
    lv1 = numeric(p)
    for(i in 1:n) lv1 = lv1 + as.numeric(object[[i]]$lv) * omega[i]
    ll0 <- sum(weight * log(lv1))
    eps=1
    iter <- 0
    scalep <- 1
    if(trace>0) cat(ll0)
    while(abs(eps)>.001 & iter<10){
        dl <- matrix(0,p,q)
        for(i in 1:n)dl <- dl + dl(object[[i]],TRUE) * omega[i]
        dl <- dl/lv1
        sc = colSums(weight * dl)
        F = crossprod(dl * weight, dl)+diag(q)*1e-6
        blub <- TRUE
        iter2 <- 0
        while(blub & iter2<10){
        thetaNew = log(theta) + scalep * solve(F, sc)
        tree$edge.length = as.numeric(exp(thetaNew))
        for(i in 1:n)object[[i]] <- update(object[[i]],tree=tree)
        lv1 = numeric(p)
        for(i in 1:n) lv1 = lv1 + as.numeric(object[[i]]$lv)  * omega[i]
        ll1 <- sum(weight * log(lv1))
        eps <- ll1 - ll0     
        if (eps < 0 || is.nan(eps)) {
            scalep = scalep/2
            eps = 1
            thetaNew = log(theta)
            ll1 = ll0
            iter2 <- iter2+1
        }
        else{
             scalep = 1;
             theta = exp(thetaNew)  
             blub=FALSE  
            }     
        }             
        iter <- iter+1
        ll0 <- ll1
    }       
    tree$edge.length <- theta
    for(i in 1:n)object[[i]] <- update(object[[i]],tree=tree)
    if(trace>0) cat("->", ll1, "\n")
    object
}


pmlMix <- function (formula, fit, m = 2, omega = rep(1/m, m), control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
{
    call <- match.call()
    form <- phangornParseFormula(formula)
    opt <- c("nni", "bf", "Q", "inv", "shape", "edge", "rate")
    optAll <- match(opt, form$left)
    optPart <- match(opt, form$right)
    AllBf <- !is.na(optAll[2])
    AllQ <- !is.na(optAll[3])
    AllInv <- !is.na(optAll[4])
    AllGamma <- !is.na(optAll[5])
    AllEdge <- !is.na(optAll[6])
    MixNni <- !is.na(optPart[1])
    MixBf <- !is.na(optPart[2])
    MixQ <- !is.na(optPart[3])
    MixInv <- !is.na(optPart[4])
    MixGamma <- !is.na(optPart[5])
    MixEdge <- !is.na(optPart[6])
    MixRate <- !is.na(optPart[7])
    if (class(fit) == "list") 
        fits <- fit
    else {
        fits <- vector("list", m) 
        for (i in 1:m) fits[[i]] <- fit
    }
    dat <- fits[[1]]$data
    p = attr(dat, "nr")
    weight = attr(dat, "weight")
    r = m
    ll = matrix(0, p, r)
    for (i in 1:r) ll[, i] = fits[[i]]$lv

    for (i in 1:r){
         pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
         fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }

    if(MixRate) rate <- rep(1,r)

    llstart = sum(weight * log(ll %*% omega))
    llold <- llstart
    ll0 <- llstart
    ll3 <- llstart
    eps0 <- 1
    iter0 <- 0
    trace = control$trace
    while (eps0 > control$eps & iter0 < control$maxit) {  #while (eps0 > 1e-6 & iter0 < 20) {
        eps1 <- 100
        iter1 <- 0
        
        if (AllQ) {
            newQ <- optimMixQ(fits, Q = fits[[1]]$Q, 
                omega = omega)[[1]]
            for (i in 1:m) fits[[i]] <- update(fits[[i]], Q = newQ)
        }
        if (AllBf) {
            newBf <- optimMixBf(fits, bf = fits[[1]]$bf, 
                omega = omega)
            for (i in 1:m) fits[[i]] <- update(fits[[i]], bf = newBf)
        }
        if (AllInv) {
            newInv <- optimMixInv(fits, inv = fits[[1]]$inv, 
                omega = omega)
            for (i in 1:m) fits[[i]] <- update(fits[[i]], Inv = newInv)
        }
        if (AllEdge) 
            fits <- optimMixEdge(fits, omega, trace=trace-1)
        for (i in 1:r) ll[, i] <- fits[[i]]$lv

        while ( abs(eps1) > 0.001 & iter1 < 3) {
             if(MixRate){
                 rate <- optimMixRate(fits, ll, weight, omega, rate)[[1]]
                 for (i in 1:r) fits[[i]] <- update(fits[[i]], rate=rate[i]) 
                 for (i in 1:r) ll[, i] <- fits[[i]]$lv
            }
            for (i in 1:r){
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
            }

            for (i in 1:r) {
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- optim.pml(fits[[i]], MixNni, MixBf, MixQ, MixInv, MixGamma, 
                    MixEdge, optRate=FALSE, control = pml.control(epsilon = 1e-8, maxit = 3,
                    trace-1), llMix = pl0, wMix = omega[i])
                 ll[, i] <- fits[[i]]$lv 

            res = optW(ll, weight, omega)
               omega = res$p
            
            if(MixRate){
                blub <- sum(rate*omega)
                rate <- rate / blub 
                tree <- fits[[1]]$tree
                tree$edge.length <-   tree$edge.length*blub
                for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, rate = rate[i])
                for (i in 1:r) ll[, i] <- fits[[i]]$lv
             }
             for (i in 1:r){
                 pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                 fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
             }
             
         }
         ll1 = sum(weight * log(ll %*% omega))
         res = optW(ll, weight, omega)
         omega = res$p
         if(MixRate){
                blub <- sum(rate*omega)
                rate <- rate / blub 
                tree <- fits[[1]]$tree
                tree$edge.length <-   tree$edge.length*blub
                for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, rate = rate[i])
                     if(trace>0) print(rate)
                     for (i in 1:r) ll[, i] <- fits[[i]]$lv
                }
         for (i in 1:r){
             pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
             fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
        }

        ll2 = sum(weight * log(ll %*% omega)) 
        eps1 = llold - ll2
        iter1 <- iter1 + 1
        llold = ll2
        }   

        ll1 <- sum(weight * log(ll %*% omega))
        eps0 <- (ll3 - ll1) / ll1
        ll3 <- ll1
        iter0 <- iter0 + 1
        if(trace>0) print(iter0)
    }
    parameter <- c(AllBf=AllBf, AllQ=AllQ, AllInv=AllInv, AllGamma=AllGamma, AllEdge=AllEdge, MixNni=MixNni, 
       MixBf=MixBf, MixQ=MixQ, MixInv=MixInv, MixGamma=MixGamma, MixEdge=MixEdge, MixRate=MixRate)
    
    df <- matrix(1, 6 ,2)
    colnames(df) <- c("#df", "group")
    rownames(df) <- c("Edge", "Shape", "Inv", "Bf", "Q", "Rate")
    df[1,1] <- length(fits[[1]]$tree$edge.length)
#    df[2,1] <- fits[[1]]$k - 1     
    df[2,1] <- fits[[1]]$k > 1
    df[3,1] <- fits[[1]]$inv > 0
    df[4,1] <- length(unique(fits[[1]]$bf)) - 1
    df[5,1] <- length(unique(fits[[1]]$Q)) - 1
    df[6,1] <- 0  
    if(MixEdge) df[1,2] = r
    if(MixGamma) df[2,2] = r
    if(MixInv) df[3,2] = r
    if(MixBf) df[4,2] = r
    if(MixQ) df[5,2] = r
    if(MixRate) df[6,1] = r-1     
    attr(logLik, "df") = sum(df[,1]*df[,2])
    converge <- c(iter=iter0, eps=eps0)
    result <- list(logLik = ll1, omega = omega, fits = fits, call = call, converge=converge, parameter=parameter, df=df)
    class(result) <- "pmlMix"
    result
}


print.pmlMix <- function(x,...){
    nc <- attr(x$fits[[1]]$data, "nc")
    nr <- attr(x$fits[[1]]$data, "nr")
    levels <- attr(x$fits[[1]]$data, "levels")
    r <- length(x$fits)   
    w <- x$fits[[1]]$weight
    w <- w[w>0] 
    type <- attr(x$fits[[1]]$data, "type")
    nc <- attr(x$fits[[1]]$data, "nc")
    ll0 = sum(w*log(w/sum(w)))

    
    bf <- matrix(0,r,nc)
    dimnames(bf) <- list(1:r, levels)
    Q <- matrix(0, r, nc*(nc-1)/2)
    dimnames(Q) <- list(1:r, NULL)

    rate <- numeric(r)
    inv <- x$fits[[1]]$inv
    shape <- numeric(r)

    for(i in 1:r){
        bf[i, ] <- x$fits[[i]]$bf
        Q[i, ] <- x$fits[[i]]$Q
        rate[i] <- x$fits[[i]]$rate
        shape[i] <- x$fits[[i]]$shape
    }
    cat("\nloglikelihood:", x$logLik, "\n")
    cat("\nunconstrained loglikelihood:", ll0, "\n") 
    cat("AIC: ", AIC(x), " BIC: ", AIC(x, k=log(nr)), "\n\n")
    cat("\nposterior:", x$omega ,"\n")   
    if(inv > 0)cat("Proportion of invariant sites:",inv,"\n")
    cat("\nRates:\n")
    cat(rate,"\n")
    cat("\nBase frequencies:  \n")
    print(bf)
    cat("\nRate matrix:\n")
    print(Q)
}


logLik.pmlMix <- function (object, ...) 
{
    res <- object$logLik
    attr(res, "df") <- sum(object$df[,1] * object$df[,2])
    class(res) <- "logLik"
    res
}
 

print.pmlPart <- function(x,...){
    df <- x$df
    nc <- attr(x$fits[[1]]$data, "nc")
    levels <- attr(x$fits[[1]]$data, "levels")
    r <- length(x$fits)   
    nc <- attr(x$fits[[1]]$data, "nc")
    nr <- attr(x$fits[[1]]$data, "nr")
    k <- x$fits[[1]]$k    

    lbf=x$df["Bf",2]
    bf <- matrix(0, lbf, nc)
    if(lbf>1)dimnames(bf) <- list(1:r, levels)
    lQ = x$df["Q",2]
    Q <- matrix(0, lQ, nc*(nc-1)/2)
    if(lQ>1)dimnames(Q) <- list(1:r, NULL)
    type <- attr(x$fits[[1]]$data, "type")
    
    loli <- numeric(r)
    rate <- numeric(r)
    shape <- numeric(r)
    sizes <- numeric(r)
    inv <- numeric(r)      
    for(i in 1:r){
        loli[i] <- x$fits[[i]]$logLik
        if(i <= lbf)bf[i, ] <- x$fits[[i]]$bf
        if(i <= lQ)Q[i, ] <- x$fits[[i]]$Q
        rate[i] <- x$fits[[i]]$rate
        shape[i] <- x$fits[[i]]$shape
        inv[i] <- x$fits[[i]]$inv
        sizes[i] <- sum(attr(x$fits[[i]]$data,"weight"))
    }
    cat("\nloglikelihood:", x$logLik, "\n")
    cat("\nloglikelihood of partitions:\n ", loli, "\n")
    cat("AIC: ", AIC(x), " BIC: ", AIC(x, k=log(sum(sizes))), "\n\n")    
    cat("Proportion of invariant sites:",inv,"\n")
    cat("\nRates:\n")
    cat(rate,"\n")
    if(k>1){
        cat("\nShape parameter:\n") 
        cat(shape,"\n")
    }
    if(type=="AA") cat("Rate matrix:",x$fits[[1]]$model, "\n")
    else{
        cat("\nBase frequencies:  \n")
        print(bf)
        cat("\nRate matrix:\n")
        print(Q)
    }
}


logLik.pmlPart <- function (object, ...) 
{
    res <- object$logLik
    attr(res, "df") <- sum(object$df[,1] * object$df[,2])
    class(res) <- "logLik"
    res
}


pmlPen <- function(object, lambda, ...){
    if(class(object)=="pmlPart") return(pmlPartPen(object, lambda,...))
    if(class(object)=="pmlMix") return(pmlMixPen(object, lambda,...))
    else stop("object has to be of class pmlPart or pmlMix")
    }
       
   
pmlPartPen <- function(object, lambda, control=pml.control(epsilon=1e-8, maxit=20, trace=1),...){
    fits <- object$fits
    
    m <- length(fits)
    K = -diag(length(fits[[1]]$tree$edge.length))
    Ktmp=K
    for(i in 1:(m-1))Ktmp = cbind(Ktmp,K)
    KM = Ktmp
    for(i in 1:(m-1))KM = rbind(KM,Ktmp)
    diag(KM) = m-1
    theta=NULL
    l = length(fits[[1]]$tree$edge.length)
    loglik = 0
    for(i in 1:m){
        theta = c(theta,fits[[i]]$tree$edge.length)
        loglik = loglik + fits[[i]]$logLik
    }
    print(loglik)
    pen = - 0.5 * lambda * t(theta)%*%KM%*%theta
    loglik = loglik - 0.5 * lambda * t(theta)%*%KM%*%theta 
    eps=1
    H  = matrix(0, m * l, m * l)
    iter=0
    trace = control$trace
    while( abs(eps)>control$eps & iter<control$maxit){
        theta=NULL
        sc = NULL
        for(i in 1:m){
            theta = c(theta,fits[[i]]$tree$edge.length)
            scoretmp = score(fits[[i]], TRUE)
            sc = c(sc,scoretmp$sc)
            H[(1:l)+l*(i-1), (1:l)+l*(i-1)] = scoretmp$F
        }
        sc = sc - lambda * KM%*% log(theta)
        thetanew = log(theta) +  solve(H + lambda*KM, sc)
        for(i in 1:m) fits[[i]]$tree$edge.length = exp(thetanew[(1:l)+(i-1)*l])
        for(i in 1:m) fits[[i]] = update.pml(fits[[i]], tree=fits[[i]]$tree)
        loglik1 = 0
        for(i in 1:m) loglik1 = loglik1 + fits[[i]]$logLik
        logLik <- loglik1
        if(trace>0)print(loglik1)
        loglik0 = loglik1
        pen = - 0.5 * lambda * t(theta)%*%KM%*%theta
        loglik1 = loglik1 - 0.5 * lambda * t(thetanew)%*%KM%*%thetanew
        eps =  (loglik - loglik1) / loglik1   
        loglik = loglik1
        theta = exp(thetanew)
        iter = iter+1
        if(trace>0)print(iter)
    }
    df = sum( diag(solve(H + lambda* KM, H)))
    
    object$df[1,1] = df
    object$df[1,2] = 1
    object$fits = fits
    object$logLik = loglik0
    attr(object$logLik, "df") = sum(object$df[,1]*object$df[,2])
    object$logLik.pen = loglik
    attr(object$logLik.pen, "df") = sum(object$df[,1]*object$df[,2])      
    object
}


pmlMixPen = function (object, lambda, optOmega=TRUE, control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
{
    fits <- object$fits
    m <- length(fits)
    K = -diag(length(fits[[1]]$tree$edge.length))
    tree <- fits[[1]]$tree
    Ktmp = K
    for (i in 1:(m - 1)) Ktmp = cbind(Ktmp, K)
    KM = Ktmp
    for (i in 1:(m - 1)) KM = rbind(KM, Ktmp)
    diag(KM) = m - 1
    theta = NULL
    l = length(fits[[1]]$tree$edge.length)
    omega <- object$omega
    dat <- fits[[1]]$data
    nr = attr(dat, "nr")
    weight = drop(attr(dat, "weight"))
    ll = matrix(0, nr, m)
    for (i in 1:m) ll[, i] = fits[[i]]$lv
    lv = drop(ll %*% omega)
    loglik = sum(weight * log(lv))
    for (i in 1:m) theta = c(theta, fits[[i]]$tree$edge.length)
    pen = - 0.5 * lambda * t(theta) %*% KM %*% theta
    loglik = loglik + pen
    print(loglik)    
    eps0 = 1 
    dl <- matrix(0, nr, m * l)
    iter0 = 0
    trace = control$trace 
    while (abs(eps0) > control$eps & iter0 < control$maxit) {
      eps = 1
      iter = 0      
      while (abs(eps) > 0.01 & iter < 5) {
        for (i in 1:m) {
            dl[, (1:l) + l * (i - 1)] <- dl(fits[[i]], TRUE) * 
                omega[i]
        }
        dl <- dl/lv
        sc = colSums(weight * dl) - lambda * KM %*% log(theta)
        H = crossprod(dl * weight, dl)
        thetanew = log(theta) + solve(H + lambda * KM, sc)
        for (i in 1:m) fits[[i]]$tree$edge.length = exp(thetanew[(1:l) + 
            (i - 1) * l])
        for (i in 1:m) {
            tree$edge.length = exp(thetanew[(1:l) + (i - 1) * l])
            fits[[i]] = update.pml(fits[[i]], tree = tree)
            ll[, i] = fits[[i]]$lv
        }
        lv = drop(ll %*% omega)
        loglik1 = sum(weight * log(lv))
        pen =  - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
        loglik1 = loglik1 + pen
        eps = abs(loglik1 - loglik)
        theta = exp(thetanew)
        loglik <- loglik1
        iter = iter + 1  
       }
       if(optOmega){
            res = optWPen(ll, weight, omega, pen)
            omega = res$p
            for (i in 1:m) {
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
                }
            } 
        lv = drop(ll %*% omega)
        loglik1 = sum(weight * log(lv))
        loglik0 =loglik1
        loglik1 = loglik1 - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
        eps0 = (loglik - loglik1) / loglik1
        theta = exp(thetanew)
        loglik <- loglik1
        iter0 = iter0 + 1
        if(trace>0) print(loglik)  
    }

    for (i in 1:m) {
        pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
        fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }
    df = sum(diag(solve(H + lambda * KM, H)))
    penalty <- list(lambda=lambda, K=KM, thetanew=thetanew, ploglik=loglik)
    object$omega = omega
    object$df[1, 1] = df
    object$df[1, 2] = 1
    object$fits = fits
    object$logLik = loglik0
    object$penalty = penalty
    object
}


optWPen = function (ll, weight, omega, pen, ...) 
{
    k = length(omega)
    nenner = 1/omega[1]
    eta = log(omega * nenner)
    eta = eta[-1]
    fn = function(eta, ll, weight, pen) {
        eta = c(0, eta)
        p = exp(eta)/sum(exp(eta))
        res = sum(weight * log(ll %*% p)) + pen
        res
    }
    if (k == 2) 
        res = optimize(f = fn, interval = c(-3, 3), lower = -3, 
            upper = 3, maximum = TRUE, tol = .Machine$double.eps^0.25, 
            ll = ll, weight = weight, pen = pen)
    else res = optim(eta, fn = fn, method = "L-BFGS-B", lower = -5, 
        upper = 5, control = list(fnscale = -1, maxit = 25), 
        gr = NULL, ll = ll, weight = weight, pen=pen)
    p = exp(c(0, res[[1]]))
    p = p/sum(p)
    result = list(par = p, value = res[[2]])
    result
} 


optNNI <- function(fit, INDEX){    
       tree = fit$tree
       ll.0 <- fit$ll.0
       loli <- fit$logLik
       bf = fit$bf
       eig = fit$eig
       k = fit$k
       w = fit$w
       g = fit$g
       rootEdges <- attr(INDEX, "root")
       .dat <- NULL
       parent = tree$edge[, 1]
       child = tree$edge[, 2]
             
       data = getCols(fit$data, tree$tip)
       datp <- rnodes(tree, data, w, g, eig, bf)       
# nicht elegant, spaeter auch raus       
       tmp = length(tree$tip.label)
       for(i in 1:length(w)).dat[i,1:tmp]=new2old.phyDat(data)       
#       datp = rnodes(fit) # raus
       
       evector <- numeric(max(parent))
       evector[child] <- tree$edge.length
       m <- dim(INDEX)[1]
       k = min(parent)
       loglik = numeric(2 * m)
       edgeMatrix <- matrix(0, 2 * m, 5)
       for (i in 1:m) {
           ei = INDEX[i, ]
           el0 = evector[INDEX[i, ]]
           l = length(datp[, 1])
           weight = fit$weight
           datn = vector("list", 4 * l)
           attr(datn, "dim") = c(l, 4)
           datn <- .dat[, ei[1:4], drop = FALSE]
           if (!(ei[5] %in% rootEdges)) 
                datn[, 1] = datp[, ei[1], drop = FALSE]
           new1 <- optim.quartet(el0[c(1, 3, 2, 4, 5)], 
               eig, bf, datn[, c(1, 3, 2, 4), drop = FALSE], g, 
               w, weight, ll.0, llcomp = fit$log)
           new2 <- optim.quartet(el0[c(1, 4, 3, 2, 5)], 
               eig, bf, datn[, c(1, 4, 3, 2), drop = FALSE], g, 
               w, weight, ll.0, llcomp = fit$log)
           loglik[(2 * i) - 1] = new1[[2]]
           loglik[(2 * i)] = new2[[2]]
           edgeMatrix[(2 * i) - 1, ] = new1[[1]]
           edgeMatrix[(2 * i), ] = new2[[1]]
           }
       list(loglik=loglik, edges = edgeMatrix)
       }


optimPartNNI <- function (object, AllEdge=TRUE,...) 
{
    tree <- object[[1]]$tree
    INDEX <- indexNNI(tree)   
    l = length(object)
    loglik0 = 0
    for(i in 1:l)loglik0 = loglik0 + logLik(object[[i]])    
    
    l = length(object)
    TMP=vector("list", l)
    for(i in 1:l){
        TMP[[i]] = optNNI(object[[i]], INDEX)
        }
    loglik=TMP[[1]][[1]] 
    for(i in 2:l)loglik=loglik+TMP[[i]][[1]]

    swap <- 0
    candidates <- loglik > loglik0

    while (any(candidates)) {
        ind = which.max(loglik)
        loglik[ind] = -Inf
        if (ind%%2) 
            swap.edge = c(2, 3)
        else swap.edge = c(2, 4)
        tree2 <- changeEdge(tree, INDEX[(ind + 1)%/%2, swap.edge], 
            INDEX[(ind + 1)%/%2, ], TMP[[1]][[2]][ind, ])
        tmpll = 0                 
        for(i in 1:l){
            if(!AllEdge)tree2 <- changeEdge(object[[i]]$tree, INDEX[(ind + 1)%/%2, swap.edge], 
                INDEX[(ind + 1)%/%2, ], TMP[[i]][[2]][ind, ]) 
            tmpll <- tmpll + update(object[[i]], tree = tree2)$logLik
            }

        if (tmpll < loglik0) 
            candidates[ind] = FALSE
        if (tmpll > loglik0) {

            swap = swap + 1
            tree <- tree2
            indi <- which(rep(colSums(apply(INDEX, 1, match, 
                INDEX[(ind + 1)%/%2, ], nomatch = 0)) > 0, each = 2))
            candidates[indi] <- FALSE
            loglik[indi] <- -Inf

            for(i in 1:l){
                if(!AllEdge)tree2 <- changeEdge(object[[i]]$tree, INDEX[(ind + 1)%/%2, swap.edge], 
                    INDEX[(ind + 1)%/%2, ], TMP[[i]][[2]][ind, ]) 
                object[[i]] <- update(object[[i]], tree = tree2)
                }
            loglik0 = 0
            for(i in 1:l)loglik0 = loglik0 + logLik(object[[i]])    
            cat(loglik0, "\n")
        }
    }
    if(AllEdge)object <- optimPartEdge(object)
    attr(object,"swap") = swap
    object
}


      
SH.test <- function (..., B = 10000, data = NULL)
{
   fits <- list(...)
   p = 1
   if (inherits(fits[[1]],"pmlPart"))
#       class(fits[[1]]) == "pmlPart") 
   {
       fits = fits[[1]]$fits
       p = length(fits)
   }
   k = length(fits)
   if (is.null(data))
       data = fits[[1]]$data
   res = NULL
   for (h in 1:p) {
       if (p > 1)
           data = fits[[h]]$data
       weight = attr(data, "weight")
       lw = length(weight)
       siteLik = matrix(0, lw, k)
       for (i in 1:k) siteLik[, i] = update(fits[[i]], data = data)$site
       ntree = k
       Lalpha <- drop(crossprod(siteLik, weight))
       Talpha <- max(Lalpha) - Lalpha
       M <- matrix(NA, k, B)
#        S <- matrix(NA, k, B)
       wvec <- rep(1L:lw, weight)
       size = length(wvec)
       for (i in 1:B) {
           boot = tabulate(sample(wvec, replace=TRUE), nbins=lw)
           M[, i] <- crossprod(siteLik, boot)
       }
       M <- M - rowMeans(M)
#        for (i in 1:B) for (j in 1:ntree) S[j, i] <- max(M[j, i] - M[, i])
       S = matrix(apply(M,2,min), k, B, byrow=TRUE)
       S = M - S
       count <- numeric(ntree)
       for (j in 1:ntree) count[j] <- sum(S[j, ] > Talpha[j])
       count <- count/B
       trees <- 1:k
       if (p == 1)
           res = cbind(trees, Lalpha, Talpha, count)
       else res = rbind(res, cbind(h, trees[-h], Lalpha[-h],
           Talpha[-h], count[-h]))
   }
   if (p == 1)
       colnames(res) <- c("Trees", "ln L", "Diff ln L", "p-value")
   else colnames(res) <- c("Partition", "Trees", "ln L", "Diff ln L",
       "p-value")
   res
}


#
# Bootstrap functions 
# multicore support
#
bootstrap.pml = function (x, bs = 100, trees = TRUE, multicore=FALSE,  ...) 
{
    data = x$data
    weight = attr(data, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for (i in 1:bs) BS[[i]] = tabulate(sample(v, replace = TRUE), 
        length(weight))
    pmlPar <- function(weights, fit, trees = TRUE, ...) {
        data = fit$data
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        fit = update(fit, data = data)
        fit = optim.pml(fit, ...)
        if (trees) {
            tree = fit$tree
            return(tree)
        }
        attr(fit, "data") = NULL
        fit
    }
    eval.success <- FALSE
    if (!eval.success & multicore) {
#  !require(parallel) ||      
        if (.Platform$GUI!="X11") {
            warning("package 'parallel' not found or GUI is used, 
            bootstrapping is performed in serial")
        } else {       
            res <- mclapply(BS, pmlPar, x, trees = trees, ...)
            eval.success <- TRUE
        } 
    }
    if (!eval.success) res <- lapply(BS, pmlPar, x, trees = trees, ...)
    if (trees) {
        class(res) = "multiPhylo"
        res = .compressTipLabel(res) # save memory
    }
    res
}


bootstrap.phyDat <- function (x, FUN, bs = 100, mc.cores=1L, ...) 
{
    weight = attr(x, "weight")
    v = rep(1:length(weight), weight)
    BS = vector("list", bs)
    for(i in 1:bs)BS[[i]]=tabulate(sample(v, replace=TRUE),length(weight)) 
    fitPar <- function(weights, data, ...){     
        ind <- which(weights > 0)
        data <- getRows(data, ind)
        attr(data, "weight") <- weights[ind]
        FUN(data,...)        
    }
    res <- mclapply(BS, fitPar, x, ..., mc.cores = mc.cores) 
    if(class(res[[1]]) == "phylo"){
        class(res) <- "multiPhylo"   
        res = .compressTipLabel(res) # save memory
    }
    res 
}


matchEdges = function(tree1, tree2){
    bp1 = bip(tree1)
    bp2 = bip(tree2)
    l = length(tree1$tip)
    fn = function(x, y){
        if(x[1]==1)return(x)
        else return(y[-x])
        } 
    bp1[] = lapply(bp1, fn, 1:l)
    bp2[] = lapply(bp2, fn, 1:l)
    match(bp1, bp2)
}


checkLabels <- function(tree, tip){
  ind <- match(tip, tree$tip.label)
  tree$tip.label <- tree$tip.label[ind]
  ind2 <- match(1:length(ind), tree$edge[, 2])
  tree$edge[ind2, 2] <- order(ind)
  tree
}


plotBS <- function (tree, BStrees, type = "unrooted", bs.col = "black", 
          bs.adj = NULL, p=80, ...) 
{
    # prop.clades raus??
    prop.clades <- function(phy, ..., part = NULL, rooted = FALSE) {
        if (is.null(part)) {
            obj <- list(...)
            if (length(obj) == 1 && class(obj[[1]]) != "phylo") 
                obj <- unlist(obj, recursive = FALSE)
            if (!identical(phy$tip, obj[[1]]$tip)) 
                obj[[1]] = checkLabels(obj[[1]], phy$tip)
            part <- prop.part(obj, check.labels = TRUE)
        }
        bp <- prop.part(phy)
        if (!rooted) {
            bp <- postprocess.prop.part(bp)
            part <- postprocess.prop.part(part)
        }
        n <- numeric(phy$Nnode)
        for (i in seq_along(bp)) {
            for (j in seq_along(part)) {
                if (identical(bp[[i]], part[[j]])) {
                    n[i] <- attr(part, "number")[j]
                    done <- TRUE
                    break
                }
            }
        }
        n
    }
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                              "unrooted", "radial"))
    if (type == "phylogram" | type == "cladogram") {
        if (!is.rooted(tree) & !is.null(tree$edge.length)) 
            tree2 = midpoint(tree)
        else tree2=tree
        plot(tree2, type = type, ...)
    }
    else plot(tree, type = type, ...)
    BStrees <- .uncompressTipLabel(BStrees)
    x = prop.clades(tree, BStrees)
    x = round((x/length(BStrees)) * 100)
    tree$node.label = x
    label = c(rep(0, length(tree$tip)), x)
    ind <- get("last_plot.phylo", envir = .PlotPhyloEnv)$edge[, 
                                                              2]
    if (type == "phylogram" | type == "cladogram") {
        root = getRoot(tree)
        label = c(rep(0, length(tree$tip)), x)
        label[root] = 0
        ind2 = matchEdges(tree2, tree)
        label = label[ind2]
        ind = which(label > p)
#        browser()
        if (is.null(bs.adj)) 
            bs.adj = c(1, 1)
        if(length(ind)>0)nodelabels(text = label[ind], node = ind, frame = "none", 
                   col = bs.col, adj = bs.adj, ...)
    }
    else {
        if (is.null(bs.adj)) 
            bs.adj = c(0.5, 0.5)
        ind2 = which(label[ind]>p)
        if(length(ind2>0))edgelabels(label[ind][ind2],ind2, frame = "none", col = bs.col, 
                   adj = bs.adj, ...)
    }
    invisible(tree)
}


pml.fit2 <- function (tree, data, bf = rep(1/length(levels), length(levels)), 
                     shape = 1, k = 1, Q = rep(1, length(levels) * (length(levels) - 1)/2), 
                     levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL, w = NULL, 
                     eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL, wMix = 0, ..., site=FALSE) 
{
    weight <- as.double(attr(data, "weight"))
    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label)) 
    k <- as.integer(k)
    m = 1
    if (is.null(eig)) 
        eig = edQt(bf = bf, Q = Q)
    if (is.null(w)) {
        w = rep(1/k, k)
        if (inv > 0) 
            w <- (1 - inv) * w
        if (wMix > 0) 
            w <- (1 - wMix) * w           
    }
    if (is.null(g)) {
        g = discrete.gamma(shape, k)
        if (inv > 0) 
            g <- g/(1 - inv)
        g <- g * rate     
    } 
#    inv0 <- inv
    if(any(g<.gEps)){
        for(i in 1:length(g)){
            if(g[i]<.gEps){
                inv <- inv+w[i]
            }
        }
        w <- w[g>.gEps]
        g <- g[g>.gEps]
#        kmax <- k
        k <- length(w)
    }
    if (is.null(INV)) 
        INV <- Matrix(lli(data, tree), sparse=TRUE)
    if (is.null(ll.0)){ 
        ll.0 <- numeric(attr(data,"nr"))    
    }
    if(inv>0)
        ll.0 <- as.matrix(INV %*% (bf * inv))              
    if (wMix > 0)
        ll.0 <- ll.0 + llMix           
    
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    root <- as.integer(node[length(node)])     
    el <- as.double(tree$edge.length)
    node = as.integer(node - nTips - 1L) #    min(node))
    edge = as.integer(edge - 1L)
    
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])    
# dlist=data, nr, nc, weight, k ausserhalb definieren  
# pmlPart einbeziehen 
    resll <- .Call("PML3", dlist=data, el, as.double(w), as.double(g), nr, nc, k, eig, as.double(bf), node, edge, nTips, root, nco, contrast, N=as.integer(length(edge))) 

    # sort(INV@i)+1L  
    ind = which(ll.0>0) # automatic in INV gespeichert

    sca = .Call("rowMax", resll, length(weight), as.integer(k)) + 1   # nr statt length(weight)
    lll = resll - sca 
    lll <- exp(lll) 
    lll <- (lll%*%w)
    lll[ind] = lll[ind] + exp(log(ll.0[ind])-sca[ind])    
    siteLik <- lll 
    siteLik <- log(siteLik) + sca
    # needs to change
    if(wMix >0) siteLik <- log(exp(siteLik) * (1-wMix) + llMix )
    loglik <- sum(weight * siteLik)
    if(!site) return(loglik)
    resll = exp(resll) 
    return(list(loglik=loglik, siteLik=siteLik, resll=resll))         
}


pml.fit4 <- function (tree, data, bf = rep(1/length(levels), length(levels)), 
                          shape = 1, k = 1, Q = rep(1, length(levels) * (length(levels) - 1)/2), 
                          levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL, w = NULL, 
                          eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL, wMix = 0, ..., site=FALSE) 
{
    weight <- as.double(attr(data, "weight"))
    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label)) 
    k <- as.integer(k)
    m = 1
    if (is.null(eig)) 
        eig = edQt(bf = bf, Q = Q)
    if (is.null(w)) {
        w = rep(1/k, k)
        if (inv > 0) 
            w <- (1 - inv) * w
        if (wMix > 0) 
            w <- (1 - wMix) * w           
    }
    if (is.null(g)) {
        g = discrete.gamma(shape, k)
        if (inv > 0) 
            g <- g/(1 - inv)
        g <- g * rate     
    } 
    #    inv0 <- inv
    if(any(g<.gEps)){
        for(i in 1:length(g)){
            if(g[i]<.gEps){
                inv <- inv+w[i]
            }
        }
        w <- w[g>.gEps]
        g <- g[g>.gEps]
        #        kmax <- k
        k <- length(w)
    }
#    .iind <- get(".iind", parent.frame())
#    .INV <- get(".INV", parent.frame())
# if(is.null(ll.0))    
    
    if (is.null(ll.0)){ 
        ll.0 <- numeric(attr(data,"nr"))    
    }
    if(inv>0)
        ll.0 <- as.matrix(INV %*% (bf * inv)) 
#    if(inv>0)  
#        ll.0 <- as.matrix(.INV %*% (bf * inv)) 
    
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    root <- as.integer(node[length(node)])     
    el <- as.double(tree$edge.length)
    node = as.integer(node - nTips - 1L) #    min(node))
    edge = as.integer(edge - 1L)
    
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])    

    siteLik <- .Call("PML4", dlist=data, el, as.double(w), as.double(g), nr, nc, k, eig, as.double(bf), 
                     node, edge, nTips, root, nco, contrast, N=as.integer(length(edge))) 
# if(inv>0) siteLik[.iind] = log(exp(siteLik[.iind]) + ll.0[.iind]) 
    ind = which(ll.0>0)
#    if(!is.null(ll.0)) siteLik[.iind] = log(exp(siteLik[.iind]) + ll.0[.iind]) 
    if(!is.null(ll.0)) siteLik[ind] = log(exp(siteLik[ind]) + ll.0[ind]) 
    if(wMix >0) siteLik <- log(exp(siteLik) * (1-wMix) + llMix )
    loglik <- sum(weight * siteLik)
    if(!site) return(loglik)
    return(list(loglik=loglik, siteLik=siteLik)) #, resll=resll         
}



pml.fit <- function (tree, data, bf = rep(1/length(levels), length(levels)), 
                     shape = 1, k = 1, Q = rep(1, length(levels) * (length(levels) - 1)/2), 
                     levels = attr(data, "levels"), inv = 0, rate = 1, g = NULL, w = NULL, 
                     eig = NULL, INV = NULL, ll.0 = NULL, llMix = NULL, wMix = 0, ..., site=FALSE) 
{
    weight <- as.double(attr(data, "weight"))
    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label)) 
    k <- as.integer(k)
    m = 1
    if (is.null(eig)) 
        eig = edQt(bf = bf, Q = Q)
    if (is.null(w)) {
        w = rep(1/k, k)
        if (inv > 0) 
            w <- (1 - inv) * w
        if (wMix > 0) 
            w <- (1 - wMix) * w           
    }
    if (is.null(g)) {
        g = discrete.gamma(shape, k)
        if (inv > 0) 
            g <- g/(1 - inv)
        g <- g * rate     
    } 
    #    inv0 <- inv
    if(any(g<.gEps)){
        for(i in 1:length(g)){
            if(g[i]<.gEps){
                inv <- inv+w[i]
            }
        }
        w <- w[g>.gEps]
        g <- g[g>.gEps]
        #        kmax <- k
        k <- length(w)
    }
    if (is.null(INV)) 
        INV <- Matrix(lli(data, tree), sparse=TRUE)
    if (is.null(ll.0)){ 
        ll.0 <- numeric(attr(data,"nr"))    
    }
    if(inv>0)
        ll.0 <- as.matrix(INV %*% (bf * inv))              
    if (wMix > 0)
        ll.0 <- ll.0 + llMix           
    
    node <- tree$edge[, 1]
    edge <- tree$edge[, 2]
    root <- as.integer(node[length(node)])     
    el <- as.double(tree$edge.length)
    node = as.integer(node - nTips - 1L) #    min(node))
    edge = as.integer(edge - 1L)
    
    contrast = attr(data, "contrast")
    nco = as.integer(dim(contrast)[1])    
    # dlist=data, nr, nc, weight, k ausserhalb definieren  
    # pmlPart einbeziehen 
    resll <- .Call("PML0", dlist=data, el, as.double(w), as.double(g), nr, nc, k, eig, as.double(bf), node, edge, nTips, root, nco, contrast, N=as.integer(length(edge))) 
    
    # sort(INV@i)+1L  
    ind = which(ll.0>0) # automatic in INV gespeichert
    
    sca = .Call("rowMax", resll, length(weight), as.integer(k)) + 1   # nr statt length(weight)
    lll = resll - sca 
    lll <- exp(lll) 
    lll <- (lll%*%w)
    lll[ind] = lll[ind] + exp(log(ll.0[ind])-sca[ind])    
    siteLik <- lll 
    siteLik <- log(siteLik) + sca
    # needs to change
    if(wMix >0) siteLik <- log(exp(siteLik) * (1-wMix) + llMix )
    loglik <- sum(weight * siteLik)
    if(!site) return(loglik)
    resll = exp(resll) 
    return(list(loglik=loglik, siteLik=siteLik, resll=resll))         
}


pml <- function (tree, data, bf = NULL, Q = NULL, inv = 0, k = 1, shape = 1, 
    rate = 1, model=NULL, ...) 
{
    call <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    pmla <- c("wMix", "llMix") 
    existing <- match(pmla, names(extras))
    wMix <- ifelse(is.na(existing[1]), 0, eval(extras[[existing[1]]], parent.frame()) )  
    llMix <- ifelse(is.na(existing[2]), 0, eval(extras[[existing[2]]], parent.frame()) )
  
    if(class(tree)!="phylo") stop("tree must be of class phylo") 
    if (is.null(attr(tree, "order")) || attr(tree, "order") == 
        "cladewise") 
        tree <- reorder(tree, "postorder")
    if(any(tree$edge.length < 0)) {
        tree$edge.length[tree$edge.length < 0] <- 1e-08
        message("negative edges length changed to 0!")
    }
    if (class(data)[1] != "phyDat") stop("data must be of class phyDat")
    if(is.null(tree$edge.length)) stop("tree must have edge weights") 
    if(any(is.na(match(tree$tip, attr(data, "names"))))) stop("tip labels are not in data")  
    data <- subset(data, tree$tip.label) # needed
    levels <- attr(data, "levels")
    weight <- attr(data, "weight")
    nr <- attr(data, "nr")
    type <- attr(data,"type")
    if(type=="AA" & !is.null(model)){
        model <- match.arg(model, .aamodels)
        getModelAA(model, bf=is.null(bf), Q=is.null(Q)) 
    }  
    if(type=="CODON") Q <- as.numeric(.syn > 0)
    if (is.null(bf)) 
        bf <- rep(1/length(levels), length(levels))
    if (is.null(Q)) 
        Q <- rep(1, length(levels) * (length(levels) - 1)/2)
    m <- 1
    eig <- edQt(bf = bf, Q = Q)

    w <- rep(1/k, k)
    if (inv > 0) 
        w <- (1 - inv) * w
    if (wMix > 0) 
        w <- wMix * w  
    g <- discrete.gamma(shape, k)
    if (inv > 0) 
        g <- g/(1 - inv)
    g <- rate * g
    inv0 <- inv
    kmax <- k    
    if(any(g<.gEps)){
        for(i in 1:length(g)){
            if(g[i]<.gEps){
                inv <- inv+w[i]
            }
        }
        w <- w[g>.gEps]
        g <- g[g>.gEps]
        k <- length(w)
    }
    
    INV <- Matrix(lli(data, tree), sparse=TRUE)
    ll.0 <- as.matrix(INV %*% (bf * inv))
    if(wMix>0) ll.0 <- ll.0 + llMix

    nr <- as.integer(attr(data, "nr")) 
    nc <- as.integer(attr(data, "nc"))
    nTips <- as.integer(length(tree$tip.label))
 
    on.exit(.C("ll_free"))
    .C("ll_init", nr, nTips, nc, as.integer(k))
    tmp <- pml.fit(tree, data, bf, shape = shape, k = k, Q = Q, 
        levels = attr(data, "levels"), inv = inv, rate = rate, g = g, w = w, 
        eig = eig, INV = INV, ll.0 = ll.0, llMix = llMix, wMix = wMix, site=TRUE) 
    
    df <- ifelse(is.ultrametric(tree), tree$Nnode, length(tree$edge.length))
    if(type=="CODON"){ 
        df <- df + (kmax>1) + (inv0 > 0) + length(unique(bf)) - 1 
        }
    else df = df + (kmax>1) + (inv0 > 0) + length(unique(bf)) - 1 + length(unique(Q)) - 1
    result = list(logLik = tmp$loglik, inv = inv, k = kmax, shape = shape,
        Q = Q, bf = bf, rate = rate, siteLik = tmp$siteLik, weight = weight, 
        g = g, w = w, eig = eig, data = data, model=model, INV = INV, 
        ll.0 = ll.0, tree = tree, lv = tmp$resll, call = call, df=df, wMix=wMix, llMix=llMix)
    if(type=="CODON"){
        result$dnds <- 1
        result$tstv <- 1
    }
    class(result) = "pml"
    result
}


optimRooted <- function(tree, data, eig=eig, w=w, g=g, bf=bf, rate=rate, ll.0=ll.0, INV=INV,
                        control = pml.control(epsilon = 1e-08, maxit = 25, trace=0), ...){
    
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08 # nicht richtig
    nTips = as.integer(length(tree$tip.label))   
    k = length(w)
    
    # optimising rooted triplets 
    optRoot0 <- function(t, tree, data, g, w, eig, bf, ll.0, k){
        l = length(tree$edge.length)
        tree$edge.length[1:(l-1)] = tree$edge.length[1:(l-1)] + t   
        tree$edge.length[l] = tree$edge.length[l] - t
        loglik = pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k) #
        loglik
    }      
    # optim edges leading to the root
    optRoot2 <- function(t, tree, data, g, w, eig, bf, ll.0, k){
        tree$edge.length = tree$edge.length + t   #c(el1+t, el2-t)
        loglik = pml.fit4(tree, data, bf=bf, g=g, w=w, eig=eig, INV=INV, ll.0=ll.0, k=k) #, INV=INV
        loglik
    }
    # scale whole tree
    scaleEdges = function(t=1, trace=0, tree, data, ...){
        fn = function(t, tree, data,...){
            tree$edge.length = tree$edge.length*t
            pml.fit4(tree, data, ...)
        }
        optimize(f=fn, interval=c(0.25,4), tree=tree, data=data, ..., maximum = TRUE,
                 tol = .00001)
    }
    parent = tree$edge[, 1]
    child = tree$edge[, 2]
    
    anc <- Ancestors(tree, 1:max(tree$edge), "parent")        
    sibs <- Siblings(tree, 1:max(tree$edge))        
    allKids <- cvector <- allChildren(tree)
    rootNode = getRoot(tree)   
    
    child2 = orderNNI(tree, nTips)   #(cvector, rootNode, nTips, TRUE)
    
    lengthList <- edgeList <- vector("list", max(tree$edge))
    for(i in tree$edge[,2]){
        pa <- anc[i]
        kids <- sibs[[i]]
        if(pa!=rootNode){
            edgeList[[i]] <- cbind(pa, c(anc[pa], kids))
            lengthList[[i]] <- c(pa, kids)
        }
        else{
            edgeList[[i]] <- cbind(pa, kids)
            lengthList[[i]] <- kids             
        }  
    }
    
    treeList <- vector("list", max(child2))
    for(i in child2){
        pa <- anc[i]
        kids <- allKids[[i]]
        treeList[[i]] <- cbind(i, c(kids, pa))
    }
    
    ll <- pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) #, INV=INV
#    if(control$trace>2)cat("ll", ll, "\n")
    eps=10
    iter = 1
    
    EL = numeric(max(tree$edge)) 
    EL[tree$edge[, 2]] = tree$edge.length
    
    eps0 =1e-8
    
    tmp <- scaleEdges(t, trace=0, tree, data, bf = bf, k=k, ll.0=ll.0, eig = eig, w=w, g=g)
#    if(control$trace>2)cat("scale", tmp[[2]], "\n")
    t = tmp[[1]]
    tree$edge.length = tree$edge.length*t        
    el = tree$edge.length
    EL[tree$edge[, 2]] = tree$edge.length
    ll2 <- pml.fit4(tree, data, bf=bf,  k=k, eig=eig, INV=INV, ll.0=ll.0, w=w, g=g) 
    
    tmptree = tree    
    
    while(eps>control$eps && iter < control$maxit){
        ll2 <- pml.fit4(tree, data, bf=bf,  k=k, eig=eig, INV=INV, ll.0=ll.0, w=w, g=g)
        loli <- rootNode
        
        children <- allKids[[rootNode]]
        kidsEl <- EL[children]  
        minEl = min(kidsEl) 
        kidsEl = kidsEl - minEl
        
        tmptree$edge = cbind(rootNode, children)
        tmptree$edge.length = kidsEl 
        
        t <- optimize(f=optRoot2,interval=c(1e-8,3), tmptree, data=data, k=k, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, maximum=TRUE)
        optRoot2(t[[1]], tmptree, data=data, k=k, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0)  
#        if(control$trace>2)cat("optRoot", t[[2]], "\n")    
        ll3 = t[[2]]
        EL[children] = kidsEl + t[[1]]     
        
        tree$edge.length = EL[tree$edge[, 2]]  
        ll2 <- pml.fit4(tree, data, bf=bf, k=k, eig=eig, INV=INV, ll.0=ll.0, w=w, g=g)  
        for(i in 1:length(child2)){ 
            dad = child2[i]
            
            #            if(dad>nTips ){ # kann raus
            pa = anc[dad]             
            while(loli != pa){                
                tmpKids= cvector[[loli]]
                tmpEdge = cbind(loli, tmpKids)
                pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf)
                loli=anc[loli] 
            }            
            pml.move(edgeList[[dad]], EL[lengthList[[dad]]], data, g, w, eig, k, nTips, bf)   
            
            children <- allKids[[dad]]
            kidsEl <- EL[children]  
            minEl = min(kidsEl) 
            maxEl = EL[dad]
            EDGE = treeList[[dad]]
            
            tmptree$edge = EDGE
            tmptree$edge.length = c(kidsEl, maxEl) 
            
            t0 = optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, k)  
            
            t = optimize(f=optRoot0, interval=c(-minEl+eps0,maxEl-eps0), tmptree, data=data, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, k=k, maximum=TRUE)
            
#            if(control$trace>2) cat("edge", t[[2]], "\n")
            if(!is.nan(t[[2]]) & t[[2]] > ll3){
                optRoot0(t[[1]], tmptree, data=data, g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, k=k)   
                EL[children] = kidsEl+t[[1]]
                EL[dad] = maxEl-t[[1]]
                ll3 = t[[2]]                   
            }
            else optRoot0(0, tmptree, data, g, w, eig, bf, ll.0, k)
            loli = dad                
            #            }
        }
        tree$edge.length = EL[tree$edge[, 2]]   
        ll2 <- pml.fit(tree, data, bf=bf, k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g)
        eps = (ll - ll2) / ll2
        
        if(control$trace>1) cat("optimRooted: ", ll, " -> ", ll2, "\n")   
        ll=ll2
        iter = iter+1
    }
    list(tree=tree, logLik=ll, c(eps=eps, iter=iter))
}



# copy node likelihoods from C to R
getNodeLogLik = function(data, i, j=1L){
    nr = attr(data, "nr")
    nc = attr(data, "nc")
    ntips = length(data)
    .Call("getLL", as.integer(i), as.integer(j-1L), as.integer(nr), as.integer(nc), as.integer(ntips))
}


# copy scaling parameters from C to R
getSC = function(data, k=1L){
    nr = attr(data, "nr")
    ntips = length(data)
    .Call("getSCM", as.integer(k),  as.integer(nr), as.integer(ntips))
}


index.nni <- function (ch, cvector, pvector, root) 
{
    p1 = pvector[ch]
    k12 = cvector[[ch]]    
    k3 = cvector[[p1]]
    k3 = k3[k3 != ch]
    kids = c(k12, k3, ch)
    parents = c(ch, ch, p1, p1)    
    if (p1 != root){    
        k4 = pvector[p1]
        kids = c(kids, k4)
        parents = c(parents, p1)
    }     
    cbind(parents, kids)
}


orderNNI <- function (tree, nTips){
    res = reorder(tree)$edge[,2]
    res = res[res>nTips]
    res
}


rooted.nni <- function(tree, data, eig, w, g, bf, rate, ll.0, INV,
    control = pml.control(epsilon = 1e-08, maxit = 25, trace=0), ...){
    ind0 = which(ll.0>0) 
    contrast = attr(data, "contrast")
    tree$edge.length[tree$edge.length < 1e-08] <- 1e-08
    nTips = as.integer(length(tree$tip.label))   
    k = length(w)
    if (is.null(attr(tree, "order")) || attr(tree, "order") == "cladewise") 
        tree <- reorder.phylo(tree, "postorder") 
    if(!is.rooted(tree))stop("tree must be rooted")
    
    attr(tree, "order") = NULL
    weight = attr(data , "weight")
    nr = as.integer(attr(data , "nr"))
    nc = as.integer(attr(data , "nc"))
    
    getEL1 <- function(t, nh){
        el = numeric(4)
        if(nh[1] > nh[2]){
            el[2] = nh[1] -nh[2]
            tnh = nh[1] + t[1]  
        }
        else{
            el[1] = nh[2] -nh[1]
            tnh = nh[2] + t[1] 
        }
        el[1:2] = el[1:2] + t[1]
        if(tnh > nh[3]) el[3] = el[3] + tnh - nh[3]
        else el[4] = el[4] - tnh + nh[3]
        el[3:4] = el[3:4] + t[2]
        el
    }
    
    optRootU <- function(t, tree, data, bf, g, w, eig, ll.0, k, INV, nh){
        tree$edge.length = getEL1(t, nh)
        pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
#        pml.fit(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
    }      

    
    getEL2 = function(t, nh){
        el = numeric(5)
        eps= 1e-6           
        nh12.min = max(nh[1:2]) + eps
        nh123.min = max(nh12.min, nh[3]) + eps
        l1 = nh[5] - nh123.min -  eps
        el[5] = l1 * t[1] + eps
        nh123 = nh[5] - el[5]
        l2 = nh123 - nh12.min - eps  
        nh12 = nh12.min + l2 * t[2]         
        el[1] = nh12 - nh[1]
        el[2] = nh12 - nh[2]
        el[3] = nh123 - nh[3]
        el[4] = nh123 - nh12
        el
    } 
    
    
    optEdgeU <- function(t, tree, data, bf, g, w, eig, ll.0, k, INV, nh){
        tree$edge.length = getEL2(t, nh)
#        pml.fit(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
        pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
    }
    
     
    child = tree$edge[, 2]   
    parent = tree$edge[, 1]
#    ll <-  pml.fit(tree, data, bf=bf, k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g)
    ll <-  pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
    llstart <- ll
    eps=.00001
    iter = 1  
    EL = numeric(max(tree$edge)) 
    EL[tree$edge[,2]] = tree$edge.length
    change = numeric(length(parent)) + 1    
    rootNode = getRoot(tree)    
    anc = Ancestors(tree, 1:max(tree$edge), "parent")  
    cvector = allChildren(tree)    
    sibs <- Siblings(tree, 1:max(tree$edge))
    
    child2 = orderNNI(tree, nTips)
    
    while(iter < 2){    
        ll2 <-  pml.fit(tree, data, bf=bf, k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g)

        nh=nodeHeight(tree)
        
        loli <- rootNode                
        pa <-rootNode    
        nchanges = 0
        ind=1
        i <- 1 
        
        tree1 <- tree2 <- tree3 <- tree
        for(i in 1:length(child2)){
            ch <- child2[i]
            dad <- anc[ch]          
            if(ch>nTips){
                
                EL[tree$edge[,2]] = tree$edge.length
                
                pa <- ifelse(dad==rootNode, rootNode ,anc[dad])   
# should avoid unnecessary movements                                                   
                while(loli != dad && loli!=rootNode){
                    if(loli==pa){ 
                        tmpKids <- sibs[[dad]]
                        tmpEdge <- cbind(pa, c(anc[pa], tmpKids))
                        pml.move(tmpEdge, EL[c(pa, tmpKids)], data, g, w, eig, k, nTips, bf)                                            
# cat("move from pa to dad \n")                                           
                        loli = dad
                    }                                    
                    else{                        
#                        cat("move loli up", loli, "dad", dad, "pa", pa, "ch", ch, "\n")
                        tmpKids = cvector[[loli]]
                        tmpEdge = cbind(loli, tmpKids)
                        pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf) 
                        loli=anc[loli]
                    }
                    
                } 
          
                if(loli == rootNode && dad!= loli){                    
                    # update all nodes
                    pml.fit(tree, data, bf=bf, k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g)                                        
#                    cat("move down loli", loli, "dad", dad, "pa", pa, "ch", ch, "\n")  
                    gd <- rev(Ancestors(tree, ch, "all")) 
                    
                    tmpKids <- sibs[[gd[2]]]
                    tmpEdge <- cbind(rootNode, tmpKids)
                    pml.move(tmpEdge, EL[tmpKids], data, g, w, eig, k, nTips, bf)                  
                    gd = gd[-1]           
                    while(length(gd)>1){
                        tmpKids <- sibs[[gd[2]]]
                        tmpEdge = cbind(gd[1], c(anc[gd[1]],tmpKids))
                        pml.move(tmpEdge, EL[c(gd[1],tmpKids)], data, g, w, eig, k, nTips, bf)      
                        gd = gd[-1]   
                    }
                    loli=dad
                }            
                
                X1 <- index.nni(ch, cvector, anc, rootNode)                 
                
#                if(loli!=rootNode){
#                    tmpKids <- c(ch, sibs[[ch]])
#                    tmpEdge <- cbind(dad, c(pa, tmpKids))
#                    tree1$edge <- tmpEdge
#                    tree1$edge.length = EL[c(dad, tmpKids)]
#                    ll0 = pml.fit(tree1, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
#                    cat("triplet", ll0, "\n")
#                }
                
                            
                if(loli!=rootNode){
                     tree1$edge <- X1
                    tree1$edge.length = abs(nh[X1[,1]] - nh[X1[,2]])
#                    ll0 = pml.fit(tree1, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
                    ll0 <- pml.fit4(tree1, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
#                    cat("quartet", ll0, ch, dad, "\n")
                }                
                
                
                if(dad == rootNode){
               
                    ll0 = pml.fit(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
                    
#                    cat("at root", ll0, ch, dad, "\n")   
                    ind2 = c(1,3,2,4)
                    ind3 = c(3,2,1,4)     
                    X2 = X3 = X1
                    X2[,2] = X1[ind2, 2]
                    X3[,2] = X1[ind3, 2]                
                    
                    tree1$edge = X1 
                    tree2$edge = X2
                    tree3$edge = X3
                    edge1 <- X1[,2]
                    edge1[4] = dad
                    res1 =optim(par = c(.1,.1), optRootU, gr=NULL, tree=tree1, data=data, nh=nh[X1[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-8, upper = 5, control = list(fnscale=-1))
                    res2 =optim(par = c(.1,.1), optRootU, gr=NULL, tree=tree2, data=data, nh=nh[X2[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-8, upper = 5, control = list(fnscale=-1))                    
                    res3 =optim(par = c(.1,.1), optRootU, gr=NULL, tree=tree3, data=data,  nh=nh[X3[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-8, upper = 5, control = list(fnscale=-1))                                       

                    ind = which.max(c(res1[[2]], res2[[2]], res3[[2]]))          
                    if(control$trace>2) cat("root", c(res1[[2]], res2[[2]], res3[[2]]), "\n")
                    
                    if(ind==1){   
                        ll2 = res1[[2]]
                        optRootU(t=res1[[1]], tree=tree1, data=data, nh=nh[X1[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k)
                        tmpEL = getEL1(res1[[1]], nh[X1[,2]])
                        tree = changeEdgeLength(tree, X1[,2], tmpEL)
                    }
                    if(ind==2){   
                        ll2 = res2[[2]]
                        optRootU(t=res2[[1]], tree=tree2, data=data, nh=nh[X2[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k) 
                        tmpEL = getEL1(res2[[1]], nh[X2[,2]])
                        tree <- changeEdge(tree, X1[c(2,3),2])                        
                        tree = changeEdgeLength(tree, X2[,2], tmpEL)
                    }
                    if(ind==2){   
                        ll2 = res3[[2]]
                        optRootU(t=res3[[1]], tree=tree3, data=data, nh=nh[X3[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k)
                        tmpEL = getEL1(res3[[1]], nh[X3[,2]])
                        tree <- changeEdge(tree, X1[c(1,3),2])
                        tree = changeEdgeLength(tree, X3[,2], tmpEL)
                    }
                }
                else{
                    loli = dad                                        
                    ind2 = c(1,3,2,4,5)
                    ind3 = c(3,2,1,4,5)
                    X2 = X3 = X1
                    X2[,2] = X1[ind2, 2]
                    X3[,2] = X1[ind3, 2]
                    tree1$edge = X1 
                    tree2$edge = X2
                    tree3$edge = X3                      
                    tt = c(.3,.5)        

                    res1 =optim(par = tt, optEdgeU, gr=NULL, tree=tree1, data, nh=nh[X1[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4, control = list(fnscale=-1))
    
                    res2 =optim(par = tt, optEdgeU, gr=NULL, tree=tree2, data, nh=nh[X2[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4, control = list(fnscale=-1))
                    
                    res3 =optim(par = tt, optEdgeU, gr=NULL, tree=tree3, data, nh=nh[X3[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k, method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4, control = list(fnscale=-1))  
                                  
                ind = which.max(c(res1[[2]], res2[[2]], res3[[2]]))     
                if(control$trace>2) cat("edge", ch, ":", c(res1[[2]], res2[[2]], res3[[2]]), "\n")    
                ll3 = max(c(res1[[2]], res2[[2]], res3[[2]]))
                
                if( (ll3 - 1e-5*ll3) < ll2){
                    loli = rootNode   
#                    ll2 <- pml.fit(tree, data, bf=bf, eig=eig, ll.0=ll.0, w=w, g=g)
                    ll2 <- pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
                    nh=nodeHeight(tree)
                    EL[tree$edge[,2]] = tree$edge.length
                    ind=0
                }   
                else{                        
                if(ind==1){   
                    ll2 = res1[[2]]
                    optEdgeU(res1[[1]], tree=tree1, data, nh=nh[X1[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k)
                    tmpEL = getEL2(res1[[1]], nh[X1[,2]])
                    tmpE = X1[,2]
                    tmpE[5] = X1[5,1]
                    tree = changeEdgeLength(tree, tmpE, tmpEL)
                }
                if(ind==2){    
                    ll2 = res2[[2]]
                    optEdgeU(res2[[1]], tree=tree2, data, nh=nh[X2[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k)
                    tmpEL = getEL2(res2[[1]], nh[X2[,2]])
                    tmpE = X2[,2]
                    tmpE[5] = X1[5,1]
                    tree <- changeEdge(tree, X1[c(2,3),2])
                    tree = changeEdgeLength(tree, tmpE, tmpEL)
                }
                if(ind==3){       
                    ll2 = res3[[2]]
                    optEdgeU(res3[[1]], tree=tree3, data, nh=nh[X3[,2]], g=g, w=w, eig=eig, bf=bf, ll.0=ll.0, INV=INV, k=k)
                    tmpEL = getEL2(res3[[1]], nh[X3[,2]])
                    tmpE = X3[,2]
                    tmpE[5] = X1[5,1]
                    tree <- changeEdge(tree, X1[c(1,3),2])
                    tree = changeEdgeLength(tree, tmpE, tmpEL)
                }
                
              }
            }
            nh=nodeHeight(tree)
            EL[tree$edge[,2]] = tree$edge.length
            loli = dad                           

            if(ind>1){
# print("NNI swap")                
                nchanges = nchanges+1 
                anc = Ancestors(tree, 1:max(tree$edge), "parent")  
                cvector = allChildren(tree)
                sibs <- Siblings(tree, 1:max(tree$edge))
                } 
            }    
            
        }
#        ll2 <- pml.fit(tree, data, bf=bf, g=g, w=w, eig=eig, ll.0=ll.0, k=k, INV=INV)
        ll2 <- pml.fit4(tree, data, bf=bf,  k=k, eig=eig, ll.0=ll.0, INV=INV, w=w, g=g) 
        
        eps = (ll - ll2) / ll2
        if(control$trace>1) cat(ll, " -> ", ll2, "\n") 
        if(control$trace>1) cat("swap:", nchanges) 
        ll=ll2
        iter = iter+1
    }
    list(tree=tree, logLik=ll, iter=iter, swap=nchanges)
}





