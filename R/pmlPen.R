#
# pmlPen penalized ML
#
pmlPen <- function(object, lambda, ...){
    if(inherits(object,"pmlPart")) return(pmlPartPen(object, lambda,...))
    if(inherits(object,"pmlMix")) return(pmlMixPen(object, lambda,...))
    else stop("object has to be of class pmlPart or pmlMix")
}


pmlPartPen <- function(object, lambda, control=pml.control(epsilon=1e-8, 
            maxit=20, trace=1),...){
    fits <- object$fits
    m <- length(fits)
    K <- -diag(length(fits[[1]]$tree$edge.length))
    Ktmp <- K
    for(i in 1:(m-1))Ktmp <- cbind(Ktmp,K)
    KM <- Ktmp
    for(i in 1:(m-1))KM <- rbind(KM,Ktmp)
    diag(KM) <- m-1
    theta <- NULL
    l <- length(fits[[1]]$tree$edge.length)
    loglik <- 0
    for(i in 1:m){
        theta <- c(theta,fits[[i]]$tree$edge.length)
        loglik <- loglik + fits[[i]]$logLik
    }
#    print(loglik)
    pen <- - 0.5 * lambda * t(theta)%*%KM%*%theta
    loglik <- loglik - 0.5 * lambda * t(theta)%*%KM%*%theta 
    eps <- 1
    H <- matrix(0, m * l, m * l)
    iter <- 0
    trace <- control$trace
    while( abs(eps)>control$eps & iter<control$maxit){
        theta <- NULL
        sc <- NULL
        for(i in 1:m){
            theta <- c(theta,fits[[i]]$tree$edge.length)
            scoretmp <- score(fits[[i]], TRUE)
            sc <- c(sc,scoretmp$sc)
            H[(1:l)+l*(i-1), (1:l)+l*(i-1)] <- scoretmp$F
        }
        sc <- sc - lambda * KM%*% log(theta)
        thetanew <- log(theta) +  solve(H + lambda*KM, sc)
        for(i in 1:m) fits[[i]]$tree$edge.length <- exp(thetanew[(1:l)+(i-1)*l])
        for(i in 1:m) fits[[i]] <- update.pml(fits[[i]], tree=fits[[i]]$tree)
        loglik1 <- 0
        for(i in 1:m) loglik1 <- loglik1 + fits[[i]]$logLik
        logLik <- loglik1
        if(trace>0)print(loglik1)
        loglik0 <- loglik1
        pen <- - 0.5 * lambda * t(theta)%*%KM%*%theta
        loglik1 <- loglik1 - 0.5 * lambda * t(thetanew)%*%KM%*%thetanew
        eps <-  (loglik - loglik1) / loglik1   
        loglik <- loglik1
        theta <- exp(thetanew)
        iter <- iter+1
        if(trace>0)print(iter)
    }
    df <- sum( diag(solve(H + lambda* KM, H)))
    
    object$df[1,1] <- df
    object$df[1,2] <- 1
    object$fits <- fits
    object$logLik <- loglik0
    attr(object$logLik, "df") <- sum(object$df[,1]*object$df[,2])
    object$logLik.pen <- loglik
    attr(object$logLik.pen, "df") <- sum(object$df[,1]*object$df[,2])      
    object
}


pmlMixPen <- function (object, lambda, optOmega=TRUE, 
                     control=pml.control(epsilon=1e-8, maxit=20, trace=1), ...) 
{
    fits <- object$fits
    m <- length(fits)
    K <- -diag(length(fits[[1]]$tree$edge.length))
    tree <- fits[[1]]$tree
    Ktmp <- K
    for (i in 1:(m - 1)) Ktmp <- cbind(Ktmp, K)
    KM <- Ktmp
    for (i in 1:(m - 1)) KM <- rbind(KM, Ktmp)
    diag(KM) <- m - 1
    theta <- NULL
    l <- length(fits[[1]]$tree$edge.length)
    omega <- object$omega
    dat <- fits[[1]]$data
    nr <- attr(dat, "nr")
    weight <- drop(attr(dat, "weight"))
    ll <- matrix(0, nr, m)
    for (i in 1:m) ll[, i] <- fits[[i]]$lv
    lv <- drop(ll %*% omega)
    loglik <- sum(weight * log(lv))
    for (i in 1:m) theta <- c(theta, fits[[i]]$tree$edge.length)
    pen <- - 0.5 * lambda * t(theta) %*% KM %*% theta
    loglik <- loglik + pen
    print(loglik)    
    eps0 <- 1 
    dl <- matrix(0, nr, m * l)
    iter0 <- 0
    trace <- control$trace 
    while (abs(eps0) > control$eps & iter0 < control$maxit) {
        eps <- 1
        iter <- 0      
        while (abs(eps) > 0.01 & iter < 5) {
            for (i in 1:m) {
                dl[, (1:l) + l * (i - 1)] <- dl(fits[[i]], TRUE) * 
                    omega[i]
            }
            dl <- dl/lv
            sc <- colSums(weight * dl) - lambda * KM %*% log(theta)
            H <- crossprod(dl * weight, dl)
            thetanew <- log(theta) + solve(H + lambda * KM, sc)
            for (i in 1:m) fits[[i]]$tree$edge.length <- exp(thetanew[(1:l) + 
                                                            (i - 1) * l])
            for (i in 1:m) {
                tree$edge.length <- exp(thetanew[(1:l) + (i - 1) * l])
                fits[[i]] <- update.pml(fits[[i]], tree = tree)
                ll[, i] <- fits[[i]]$lv
            }
            lv <- drop(ll %*% omega)
            loglik1 <- sum(weight * log(lv))
            pen <-  - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
            loglik1 <- loglik1 + pen
            eps <- abs(loglik1 - loglik)
            theta <- exp(thetanew)
            loglik <- loglik1
            iter <- iter + 1  
        }
        if(optOmega){
            res <- optWPen(ll, weight, omega, pen)
            omega <- res$p
            for (i in 1:m) {
                pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
                fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
            }
        } 
        lv <- drop(ll %*% omega)
        loglik1 <- sum(weight * log(lv))
        loglik0 <- loglik1
        loglik1 <- loglik1 - 0.5 * lambda * t(thetanew) %*% KM %*% thetanew
        eps0 <- (loglik - loglik1) / loglik1
        theta <- exp(thetanew)
        loglik <- loglik1
        iter0 <- iter0 + 1
        if(trace>0) print(loglik)  
    }
    
    for (i in 1:m) {
        pl0 <- ll[, -i, drop = FALSE] %*% omega[-i]
        fits[[i]] <- update(fits[[i]], llMix = pl0, wMix = omega[i])
    }
    df <- sum(diag(solve(H + lambda * KM, H)))
    penalty <- list(lambda=lambda, K=KM, thetanew=thetanew, ploglik=loglik)
    object$omega <- omega
    object$df[1, 1] <- df
    object$df[1, 2] <- 1
    object$fits <- fits
    object$logLik <- loglik0
    object$penalty <- penalty
    object
}


optWPen <- function (ll, weight, omega, pen, ...) 
{
    k <- length(omega)
    nenner <- 1/omega[1]
    eta <- log(omega * nenner)
    eta <- eta[-1]
    fn <- function(eta, ll, weight, pen) {
        eta <- c(0, eta)
        p <- exp(eta)/sum(exp(eta))
        res <- sum(weight * log(ll %*% p)) + pen
        res
    }
    if (k == 2) 
        res <- optimize(f = fn, interval = c(-3, 3), lower = -3, upper = 3, 
                maximum = TRUE, tol = .Machine$double.eps^0.25, 
                ll = ll, weight = weight, pen = pen)
    else res <- optim(eta, fn = fn, method = "L-BFGS-B", lower = -5, 
                     upper = 5, control = list(fnscale = -1, maxit = 25), 
                     gr = NULL, ll = ll, weight = weight, pen=pen)
    p <- exp(c(0, res[[1]]))
    p <- p/sum(p)
    result <- list(par = p, value = res[[2]])
    result
} 
