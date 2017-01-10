
compatible_2 <- function(obj1, obj2) 
{
    ntaxa <- length(obj1$labels)
    msk <- !as.raw(2^(8 - (ntaxa%%8)) - 1)
    r0 <- as.raw(0)
    arecompatible2 <- function (x, y, msk, r0) {
        foo <- function(v) {
            lv <- length(v)
            v[lv] <- v[lv] & msk
            as.integer(all(v == r0))
        }
        nE <- foo(x & y) + foo(x & !y) + foo(!x & y) + foo(!x & !y)
        if (nE > 0) TRUE
        else FALSE
    }  
    m1 <- obj1$matsplit
    m2 <- obj2$matsplit
    n1 <- ncol(m1)
    n2 <- ncol(m2)
    res <- rep(TRUE, n1)
    for (i in 1:n1){
        j=1
        while(j <= n2){
            if (!arecompatible2(m1[, i], m2[, j], msk, r0)){
                res[i]=FALSE
                break()
            }
            j=j+1L
        } 
    }    
    res
}          


addConfidences_MultiPhylo <- function(spl, trees){
    
    fun <- function(spl, intersect_labels){
        spl2 <- spl
        
        index <- match(attr(spl, "labels"), intersect_labels)
        attr(spl2, "labels") <- intersect_labels 
        for(i in 1:length(spl2)){
            spl2[[i]] <- sort(na.omit(index[spl[[i]]]))
        }
        l_spl <- lengths(spl2)
        l <- length(intersect_labels)
        ind <- which( (l_spl > 1) & (l_spl < (l-1L)) ) 
        if(length(ind)==0)return(NULL)
        list(spl=spl2[ind], index=ind)
    }
    
#    spl_net <- net$splits
    spl_labels <- attr(spl, "labels")
    zaehler <- numeric(length(spl))
    nenner <- numeric(length(spl))
    
    for(i in 1:length(trees)){
        print(i)
        intersect_labels <- intersect(trees[[i]]$tip.label, spl_labels)
        if(length(intersect_labels) > 3){    
            tmp <- fun(spl, intersect_labels)
            if(!is.null(tmp)){
                tree_spl <- as.splits(trees[[i]])
                if(!identical(intersect_labels, trees[[i]]$tip.label))
                    tree_spl <- fun(tree_spl, intersect_labels)[[1]]
                comp <- compatible_2(as.bitsplits(tmp[[1]]), as.bitsplits(tree_spl))
                #            print(comp)
                ind <- tmp$index
                zaehler[ind] <- zaehler[ind] + comp 
                nenner[ind] <- nenner[ind] + 1L
            }
        }    
        #        print(zaehler)
    }
    confidences <- zaehler / nenner
    attr(spl, "confidences") <- confidences
    spl
}            

