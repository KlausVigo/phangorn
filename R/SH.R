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