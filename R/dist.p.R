
dist.p <- function (x, cost="polymorphism", ignore.indels=TRUE) 
{
    if (class(x) != "phyDat") 
        stop("x has to be element of class phyDat")

    l = length(x)
    weight <- attr(x, "weight")
    n <- length(attr(x, "allLevels"))
    d = numeric((l * (l - 1))/2)
    lev = attr(x, "allLevels")    
    if(is.null(cost)){ 
        cost <- 1 - diag(n)
        dimnames(cost) = list(lev, lev)
    }    
#    if(cost=="polymorphism" && attr(x, "type")=="DNA"){   
    if(cost=="polymorphism"){
        costLev = c('a','c','t','u','g','x','m','r','w','s','y','k','v','h','d','b','-','?','n')
        
    cost <- matrix(c(
       #a,c,t,u,g,X,m,r,w,s,y,k,v,h,d,b,-,?,n,
        0,2,2,2,2,1,1,1,1,3,3,3,2,2,2,4,2,0,0, #a
        2,0,2,2,2,1,1,3,3,1,1,3,2,2,4,2,2,0,0, #c
        2,2,0,0,2,1,3,3,1,3,1,1,4,2,2,2,2,0,0, #t
        2,2,0,0,2,1,3,3,1,3,1,1,4,2,2,2,2,0,0, #u
        2,2,2,2,0,1,3,1,3,1,3,1,2,4,2,2,2,0,0, #g
        1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0, #X
        1,1,3,3,3,1,0,2,2,2,2,4,1,1,3,3,3,0,0, #m
        1,3,3,3,1,1,2,0,2,2,4,2,1,3,1,3,3,0,0, #r
        1,3,1,1,3,1,2,2,0,4,2,2,3,1,1,3,3,0,0, #w
        3,1,3,3,1,1,2,2,4,0,2,2,1,3,3,1,3,0,0, #s
        3,1,1,1,3,1,2,4,2,2,0,2,3,1,3,1,3,0,0, #y
        3,3,1,1,1,1,4,2,2,2,2,0,3,3,1,1,3,0,0, #k
        2,2,4,4,2,1,1,1,3,1,3,3,0,2,2,2,4,0,0, #v
        2,2,2,2,4,1,1,3,1,3,1,3,2,0,2,2,4,0,0, #h
        2,4,2,2,2,1,3,1,1,3,3,1,2,2,0,2,4,0,0, #d
        4,2,2,2,2,1,3,3,3,1,1,1,2,2,2,0,4,0,0, #b
        2,2,2,2,2,1,3,3,3,3,3,3,4,4,4,4,0,0,0, #-
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #?
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),#n
           ncol = 19,nrow=19,dimnames=list(costLev,costLev))
    }
    
    lev1 = dimnames(cost)[[1]]
    

    if(any(is.na(match(lev, lev1)))) stop("Levels of x are not in levels of cost matrix!")

        if (ignore.indels) {
            cost["-",]=0
            cost[,"-"]=0
        } 

    
    cost <- cost[lev, lev]
    
        k = 1
        for (i in 1:(l - 1)) {
            for (j in (i + 1):l) {
                d[k] = sum(weight * cost[cbind(x[[i]], x[[j]])])
                k = k + 1
            }
        }
    attr(d, "Size") <- l
    if (is.list(x)) 
        attr(d, "Labels") <- names(x)
    else attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- "p"
    class(d) <- "dist"
    return(d)
}





