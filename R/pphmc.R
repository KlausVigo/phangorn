# Probabilistic Path Hamiltonian Monte Carlo 

# translated from Logprior.py

phyloLogprior_exp <- function(tree, scale=0.1, grad=FALSE){
    edge <- tree$edge.length
    if( !grad ) return( - sum(edge)/scale - log(scale) * length(edge) )
    rep(-1, length(edge))/scale 
}


# translated from Logposterior.py

# a mollifier for surrogate construction
mollifier <- function(x, delta){
    x[x<delta] <- 0.5 / delta *(x[x<delta]*x[x<delta]+delta*delta)
    x
}


# negative log-posterior
Logpost <- function(tree, data, bf, eig, scale=0.1, surrogate=FALSE, delta = 0.01, ...){
    if(surrogate)
        tree$edge.length <- mollifier(tree$edge.length, delta) 
    -pml.fit(tree, data, bf=bf, eig=eig, ...) - phyloLogprior_exp(tree, scale)
}


# tree, branch, D, U, U_inv, pden, L, scale=0.1,

#GradLogpost <- function(tree, data, bf, eig, scale=0.1, surrogate=FALSE, delta = 0.01){
#    if(surrogate){
#        tree$edge.length <- mollifier(tree$edge.length, delta) 
#        return (-phyloLoglikelihood(tree, maped_branch, D, U, U_inv, pden, L, grad=True)  
#        - phyloLogprior_exp(tree, scale, grad=TRUE)) * np.minimum(branch,delta)/delta
#    }
#    -phyloLoglikelihood(tree, branch, D, U, U_inv, pden, L, grad=True) - phyloLogprior_exp(tree, scale, grad=TRUE)
#}    


# score function / gradient?
score_2 <- function (fit, transform=FALSE) 
{
    tree <- fit$tree
    child <- tree$edge[, 2]
    l <- length(child)
    sc <- numeric(l)
    weight <- as.numeric(fit$weight)
    f <- drop(exp(fit$site))
    dl <- dl(fit, transform)
    dl <- dl/f
    sc <- colSums(weight * dl)
    sc
}



    