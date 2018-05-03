#
# pmlMix
#
optimMixQ <- function(object, Q=c(1, 1, 1, 1, 1, 1), omega,...){
    l <- length(Q)
    Q <- Q[-l]
    Q <- sqrt(Q)
    fn <- function(Q, object, omega,...) {
        Q <- c(Q^2, 1)
        weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + 
            as.numeric(update(object[[i]], Q=Q, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res <- optim(par=Q, fn=fn, gr=NULL, method="L-BFGS-B", lower=0, 
                upper=Inf, control=list(fnscale = -1, maxit=25), 
                object=object, omega=omega,...)
    res[[1]] <- c(res[[1]]^2, 1)
    res
}


optimMixBf <- function(object, bf=c(.25,.25,.25,.25), omega,...){
    l <- length(bf)
    nenner <- 1/bf[l]
    lbf <- log(bf * nenner)
    lbf <- lbf[-l]
    fn <- function(lbf, object, omega,...) {
        bf <- exp(c(lbf,0))
        bf <- bf/sum(bf)
        weight <- object[[1]]$weight
        n <- length(omega)
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n)result <- result + 
            as.numeric(update(object[[i]], bf=bf, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res <- optim(par=lbf, fn=fn, gr=NULL, method="Nelder-Mead", 
                control=list(fnscale=-1, maxit=500), object, omega=omega,...)
    #    print(res[[2]])
    bf <- exp(c(res[[1]],0))
    bf <- bf/sum(bf)
}


optimMixInv <- function(object, inv=0.01, omega,...){
    fn <- function(inv, object, omega,...) {
        n <- length(omega)
        weight <- object[[1]]$weight
        p <- length(weight)
        result <- numeric(p)
        for(i in 1:n) result <- result + 
            as.numeric(update(object, inv=inv, ...)$lv) * omega[i]
        result <- sum(weight %*% log(result))
        result 
    }
    res <- optimize(f=fn, interval=c(0,1), lower=0, upper=1, maximum=TRUE,
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
    ui <- diag(r-1)
    ui <- rbind(-omega[-r], ui)
    ci <- c(-1, rep(0, r-1))
    res <- constrOptim(rate0, fn, grad=NULL, ui=ui, ci=ci, mu=1e-04, 
            control=list(fnscale=-1), method="Nelder-Mead", outer.iterations=50,
            outer.eps = 1e-05, fits=fits, ll=ll, weight=weight, omega=omega)
    rate <- res[[1]]
    res[[1]] <- c(rate, (1- sum(rate *omega[-r]))/omega[r])
    res
}


optW <- function (ll, weight, omega,...) 
{
    k <- length(omega)
    nenner <- 1/omega[1]
    eta <- log(omega * nenner)
    eta <- eta[-1]
    fn <- function(eta, ll, weight) {
        eta <- c(0,eta)
        p <- exp(eta)/sum(exp(eta))
        res <- sum(weight * log(ll %*% p)) 
        res
    }
    if(k==2)res <- optimize(f=fn , interval=c(-3,3) , lower=-3, upper=3, 
               maximum=TRUE, tol=.Machine$double.eps^0.25, ll=ll, weight=weight) 
    else res <- optim(eta, fn=fn, method="L-BFGS-B", lower=-5, upper=5, 
            control = list(fnscale=-1, maxit=25), gr=NULL, ll=ll, weight=weight)
    
    p <- exp(c(0,res[[1]]))
    p <- p/sum(p)
    result <- list(par = p, value = res[[2]])
    result
}


optimMixEdge <- function(object, omega, trace=1,...){
    tree <- object[[1]]$tree
    theta <- object[[1]]$tree$edge.length
    weight <- as.numeric(attr(object[[1]]$data,"weight"))
    n <- length(omega)
    p <- length(weight)
    q <- length(theta)
    lv1 <- numeric(p)
    for(i in 1:n) lv1 <- lv1 + as.numeric(object[[i]]$lv) * omega[i]
    ll0 <- sum(weight * log(lv1))
    eps <- 1
    iter <- 0
    scalep <- 1
    if(trace>0) cat(ll0)
    while(abs(eps)>.001 & iter<10){
        dl <- matrix(0,p,q)
        for(i in 1:n)dl <- dl + dl(object[[i]],TRUE) * omega[i]
        dl <- dl/lv1
        sc <- colSums(weight * dl)
        F <- crossprod(dl * weight, dl)+diag(q)*1e-6
        blub <- TRUE
        iter2 <- 0
        while(blub & iter2<10){
            thetaNew <- log(theta) + scalep * solve(F, sc)
            tree$edge.length <- as.numeric(exp(thetaNew))
            for(i in 1:n)object[[i]] <- update(object[[i]],tree=tree)
            lv1 <- numeric(p)
            for(i in 1:n) lv1 <- lv1 + as.numeric(object[[i]]$lv)  * omega[i]
            ll1 <- sum(weight * log(lv1))
            eps <- ll1 - ll0     
            if (eps < 0 || is.nan(eps)) {
                scalep <- scalep/2
                eps <- 1
                thetaNew <- log(theta)
                ll1 <- ll0
                iter2 <- iter2+1
            }
            else{
                scalep <- 1
                theta <- exp(thetaNew)  
                blub <- FALSE  
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




#' Phylogenetic mixture model
#' 
#' Phylogenetic mixture model.
#' 
#' The \code{formula} object allows to specify which parameter get optimized.
#' The formula is generally of the form \code{edge + bf + Q ~ rate + shape +
#' \dots{}}, on the left side are the parameters which get optimized over all
#' mixtures, on the right the parameter which are optimized specific to each
#' mixture. The parameters available are \code{"nni", "bf", "Q", "inv",
#' "shape", "edge", "rate"}.  Each parameters can be used only once in the
#' formula.  \code{"rate"} and \code{"nni"} are only available for the right
#' side of the formula. On the other hand parameters for invariable sites are
#' only allowed on the left-hand side.  The convergence of the algorithm is
#' very slow and is likely that the algorithm can get stuck in local optima.
#' 
#' @aliases pmlMix pmlPen
#' @param formula a formula object (see details).
#' @param fit an object of class \code{pml}.
#' @param m number of mixtures.
#' @param omega mixing weights.
#' @param control A list of parameters for controlling the fitting process.
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{pmlMix} returns a list with elements
#' \item{logLik}{log-likelihood of the fit} \item{omega}{mixing weights.}
#' \item{fits}{fits for the final mixtures.}
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{pml}},\code{\link{pmlPart}},\code{\link{pmlCluster}}
#' @keywords cluster
#' @examples
#' 
#' \dontrun{
#' X <- allSitePattern(5)
#' tree <- read.tree(text = "((t1:0.3,t2:0.3):0.1,(t3:0.3,t4:0.3):0.1,t5:0.5);")
#' fit <- pml(tree,X, k=4)
#' weights <- 1000*exp(fit$site) 
#' attr(X, "weight") <- weights
#' fit1 <- update(fit, data=X, k=1)
#' fit2 <- update(fit, data=X)
#' 
#' (fitMixture <- pmlMix(edge~rate, fit1 , m=4))
#' (fit2 <- optim.pml(fit2, optGamma=TRUE))
#' 
#' 
#' data(Laurasiatherian)
#' dm <- dist.logDet(Laurasiatherian)
#' tree <- NJ(dm)
#' fit <- pml(tree, Laurasiatherian)
#' fit <- optim.pml(fit)
#' 
#' fit2 <- update(fit, k=4)
#' fit2 <- optim.pml(fit2, optGamma=TRUE)
#' 
#' fitMix <- pmlMix(edge ~ rate, fit, m=4)
#' fitMix
#' 
#' 
#' #
#' # simulation of mixture models 
#' #
#' \dontrun{
#' X <- allSitePattern(5)
#' tree1 <- read.tree(text = "((t1:0.1,t2:0.5):0.1,(t3:0.1,t4:0.5):0.1,t5:0.5);")
#' tree2 <- read.tree(text = "((t1:0.5,t2:0.1):0.1,(t3:0.5,t4:0.1):0.1,t5:0.5);")
#' tree1 <- unroot(tree1)
#' tree2 <- unroot(tree2)
#' fit1 <- pml(tree1,X)
#' fit2 <- pml(tree2,X)
#' 
#' weights <- 2000*exp(fit1$site) + 1000*exp(fit2$site)
#' attr(X, "weight") <- weights
#' 
#' fit1 <- pml(tree1, X)
#' fit2 <- optim.pml(fit1)
#' logLik(fit2)
#' AIC(fit2, k=log(3000))
#' 
#' fitMixEdge <- pmlMix( ~ edge, fit1, m=2)
#' logLik(fitMixEdge)
#' AIC(fitMixEdge, k=log(3000))
#' 
#' fit.p <- pmlPen(fitMixEdge, .25)
#' logLik(fit.p)
#' AIC(fit.p, k=log(3000))
#' }
#' }
#' 
#' @export pmlMix
pmlMix <- function (formula, fit, m = 2, omega = rep(1/m, m), 
                    control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
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
    if (inherits(fit,"pmlMix")) fits <- fit$fits     
    if (inherits(fit,"list")) 
        fits <- fit
    if (inherits(fit,"pml")){
        fits <- vector("list", m) 
        for (i in 1:m) fits[[i]] <- fit
    }
    dat <- fits[[1]]$data
    p <- attr(dat, "nr")
    weight <- attr(dat, "weight")
    r <- m
    ll <- matrix(0, p, r)
    for (i in 1:r) ll[, i] <- fits[[i]]$lv
    
    for (i in 1:r){
        pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
        fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }
    
    if(MixRate) rate <- rep(1,r)
    
    llstart <- sum(weight * log(ll %*% omega))
    llold <- llstart
    ll0 <- llstart
    ll3 <- llstart
    eps0 <- 1
    iter0 <- 0
    trace <- control$trace
    while (eps0 > control$eps & iter0 < control$maxit) {  
        #while (eps0 > 1e-6 & iter0 < 20) {
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
                fits[[i]] <- optim.pml(fits[[i]], MixNni, MixBf, MixQ, MixInv, 
                    MixGamma, MixEdge, optRate=FALSE, control = 
                    pml.control(epsilon = 1e-8, maxit = 3, trace-1), 
                    llMix = pl0, wMix = omega[i])
                ll[, i] <- fits[[i]]$lv 
                
                res <- optW(ll, weight, omega)
                omega <- res$p
                
                if(MixRate){
                    blub <- sum(rate*omega)
                    rate <- rate / blub 
                    tree <- fits[[1]]$tree
                    tree$edge.length <-   tree$edge.length*blub
                    for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, 
                                rate = rate[i])
                    for (i in 1:r) ll[, i] <- fits[[i]]$lv
                }
                for (i in 1:r){
                    pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                    fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
                }
                
            }
            ll1 <- sum(weight * log(ll %*% omega))
            res <- optW(ll, weight, omega)
            omega <- res$p
            if(MixRate){
                blub <- sum(rate*omega)
                rate <- rate / blub 
                tree <- fits[[1]]$tree
                tree$edge.length <-   tree$edge.length*blub
                for (i in 1:r) fits[[i]]<-update(fits[[i]], tree=tree, 
                                                 rate = rate[i])
                if(trace>0) print(rate)
                for (i in 1:r) ll[, i] <- fits[[i]]$lv
            }
            for (i in 1:r){
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
            }
            
            ll2 <- sum(weight * log(ll %*% omega)) 
            eps1 <- llold - ll2
            iter1 <- iter1 + 1
            llold <- ll2
        }   
        
        ll1 <- sum(weight * log(ll %*% omega))
        eps0 <- (ll3 - ll1) / ll1
        ll3 <- ll1
        iter0 <- iter0 + 1
        if(trace>0) print(iter0)
    }
    parameter <- c(AllBf=AllBf, AllQ=AllQ, AllInv=AllInv, AllGamma=AllGamma, 
          AllEdge=AllEdge, MixNni=MixNni, MixBf=MixBf, MixQ=MixQ, MixInv=MixInv,
          MixGamma=MixGamma, MixEdge=MixEdge, MixRate=MixRate)
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
    if(MixEdge) df[1,2] <- r
    if(MixGamma) df[2,2] <- r
    if(MixInv) df[3,2] <- r
    if(MixBf) df[4,2] <- r
    if(MixQ) df[5,2] <- r
    if(MixRate) df[6,1] <- r-1     
    attr(logLik, "df") <- sum(df[,1]*df[,2])
    converge <- c(iter=iter0, eps=eps0)
    result <- list(logLik = ll1, omega = omega, fits = fits, call = call, 
                   converge=converge, parameter=parameter, df=df)
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
    ll0 <- sum(w*log(w/sum(w)))
    
    
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

