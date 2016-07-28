
#
# add codon models, change to phyDat statt 3* 
#
simSeq <- function (x, ...) 
    UseMethod("simSeq")


simSeq.phylo = function(x, l=1000, Q=NULL, bf=NULL, rootseq=NULL, type = "DNA", model="USER",
                  levels = NULL, rate=1, ancestral=FALSE, ...){
    
    pt <- match.arg(type, c("DNA", "AA", "USER", "CODON"))
    if (pt == "DNA") 
        levels <- c("a", "c", "g", "t")
    if (pt == "AA") 
        levels <- c("a", "r", "n", "d", "c", "q", "e", "g", "h", "i", 
                    "l", "k", "m", "f", "p", "s", "t", "w", "y", "v")
    if (pt == "CODON"){
        levels <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act", 
          "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att", 
          "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", "cga", 
          "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac", 
          "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg", 
          "ggt", "gta", "gtc", "gtg", "gtt", "tac", "tat", 
          "tca", "tcc", "tcg", "tct", "tgc", "tgg", "tgt", "tta", 
          "ttc", "ttg", "ttt")
        Q <- as.numeric(.syn > 0)
    }
    if (pt == "USER") 
        if(is.null(levels))stop("levels have to be supplied if type is USER")
    
    lbf = length(levels)
    
    if (type == "AA" & !is.null(model)) {
        #        model <- match.arg(model, c("USER", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        model <- match.arg(model, c("USER", .aamodels))
        if(model!="USER")getModelAA(model, bf=is.null(bf), Q=is.null(Q))
    }
    
    if(is.null(bf)) bf = rep(1/lbf,lbf)
    if(is.null(Q)) Q = rep(1,lbf*(lbf-1)/2)
    if(is.matrix(Q)) Q=Q[lower.tri(Q)]
    eig = edQt(Q, bf)
    
    m = length(levels)    
    
    if(is.null(rootseq))rootseq = sample(levels, l, replace=TRUE, prob=bf)
    x = reorder(x) 
    edge = x$edge
    nNodes = max(edge)
    res = matrix(NA, l, nNodes)
    parent <- as.integer(edge[, 1])
    child <- as.integer(edge[, 2])
    root <- as.integer(parent[!match(parent, child, 0)][1])  
    res[, root] = rootseq   
    tl = x$edge.length
    for(i in 1:length(tl)){
        from = parent[i] 
        to = child[i]
        P = getP(tl[i], eig, rate)[[1]]
        for(j in 1:m){
            ind = res[,from]==levels[j]
            res[ind,to] = sample(levels, sum(ind), replace=TRUE, prob=P[,j])
        }
    }
    k = length(x$tip.label)
    label = c(x$tip.label, as.character((k+1):nNodes))
    colnames(res)=label 
    if(!ancestral)res = res[, x$tip.label, drop=FALSE]
    if(pt=="DNA") return(phyDat.DNA(as.data.frame(res, stringsAsFactors = FALSE), return.index=TRUE))
    if(pt=="AA") return(phyDat.AA(as.data.frame(res, stringsAsFactors = FALSE), return.index=TRUE))
    if(pt=="USER") return(phyDat.default(as.data.frame(res, stringsAsFactors = FALSE), levels = levels, return.index=TRUE))
    if(pt=="CODON"){ 
        res <- apply(res, 2, function(x)unlist(strsplit(x, "")))
        return(phyDat.codon(as.data.frame(res, stringsAsFactors = FALSE)))
    }
}        



simSeq.pml <- function(x, ancestral=FALSE, ...){
    g = x$g
    w = x$w
    if(x$inv>0){
        w = c(x$inv, w)
        g = c(0.0, g)
    }
    n = length(w)
    res = vector("list", n)
    y = sample(n, sum(x$weight), replace=TRUE, prob=w)
    levels = attr(x$data, "levels")
    type = attr(x$data, "type")
    for(i in 1:n){
        l = sum(y==i)
        res[[i]] = simSeq(x$tree, l, Q=x$Q, bf=x$bf, type=type, levels=levels, rate=g[i], ancestral=ancestral)  
    }
    x = call("c.phyDat", quote(res[[1]]))
    if(n>1) x <- parse(text= paste("c(", "res[[1]]", paste0(",res[[", 2:n, "]]", collapse=""), ")"))
    eval(x)    
}
