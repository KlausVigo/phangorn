#
# add codon models, change to phyDat statt 3* 
#


#' Simulate sequences.
#' 
#' Simulate sequences for a given evolutionary tree.
#' 
#' \code{simSeq} is now a generic function to simulate sequence alignments.  It
#' is quite flexible and allows to generate DNA, RNA, amino acids or binary
#' sequences.  It is possible to give a \code{pml} object as input simSeq
#' return a \code{phyDat} from these model.  There is also a more low level
#' version, which lacks rate variation, but one can combine different
#' alignments having their own rate (see example). The rate parameter acts like
#' a scaler for the edge lengths.
#' 
#' @param x a phylogenetic tree \code{tree}, i.e. an object of class
#' \code{phylo} or and object of class \code{pml}.
#' @param l length of the sequence to simulate.
#' @param Q the rate matrix.
#' @param bf base frequencies.
#' @param rootseq a vector of length l containing the root sequence, other root
#' sequence is randomly generated.
#' @param type Type of sequences ("DNA", "AA" or "USER").
#' @param model Amino acid models: e.g. "WAG", "JTT", "Dayhoff" or "LG"
#' @param levels \code{levels} takes a character vector of the different bases,
#' default is for nucleotide sequences, only used when type = "USER".
#' @param rate mutation rate or scaler for the edge length, a numerical value
#' greater than zero.
#' @param ancestral Return ancestral sequences?
#' @param \dots Further arguments passed to or from other methods.
#' @return \code{simSeq} returns an object of class phyDat.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{phyDat}}, \code{\link{pml}}, \code{\link{SOWH.test}}
#' @keywords cluster
#' @examples
#' 
#' \dontrun{
#' data(Laurasiatherian)
#' tree = nj(dist.ml(Laurasiatherian))
#' fit = pml(tree, Laurasiatherian, k=4)
#' fit = optim.pml(fit, optNni=TRUE, model="GTR", optGamma=TRUE)
#' data = simSeq(fit)
#' }
#' 
#' tree = rtree(5)
#' plot(tree)
#' nodelabels()
#' 
#' # Example for simple DNA alignment
#' data = simSeq(tree, l = 10, type="DNA", bf=c(.1,.2,.3,.4), Q=1:6)
#' as.character(data)
#' 
#' # Example to simulate discrete Gamma rate variation
#' rates = discrete.gamma(1,4)
#' data1 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[1])
#' data2 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[2])
#' data3 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[3])
#' data4 = simSeq(tree, l = 100, type="AA", model="WAG", rate=rates[4])
#' data <- c(data1,data2, data3, data4)
#' 
#' write.phyDat(data, file="temp.dat", format="sequential",nbcol = -1, colsep = "")
#' unlink("temp.dat") 
#' 
#' @rdname simSeq
#' @export simSeq
simSeq <- function (x, ...) 
    UseMethod("simSeq")


#' @rdname simSeq
#' @method simSeq phylo
#' @export
simSeq.phylo = function(x, l=1000, Q=NULL, bf=NULL, rootseq=NULL, type = "DNA", model=NULL,
                  levels = NULL, rate=1, ancestral=FALSE, ...){
   
    
    if (!is.null(model)) {
        #        model <- match.arg(model, c("USER", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24"))
        model <- match.arg(model, .aamodels) #match.arg(model, c("USER", .aamodels)) 
        getModelAA(model, bf=is.null(bf), Q=is.null(Q))
        type = "AA"
    }
     
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
    for(i in seq_along(tl)){
        from = parent[i] 
        to = child[i]
        P = getP(tl[i], eig, rate)[[1]]
        # avoid numerical problems for larger P and small t        
        if(any(P < 0)) P[P<0] = 0
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


#' @rdname simSeq
#' @method simSeq pml
#' @export
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
