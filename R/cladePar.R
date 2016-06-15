cladePar = function(tree, node, edge.color="red", tip.color=edge.color, edge.width = 1, edge.lty = 1, x=NULL, plot=FALSE, ...){
    if(is.null(x)){
        m = max(tree$edge)
        x=list(edge=data.frame(color=rep("black",m), width = rep(1, m), lty =  rep(1, m), stringsAsFactors = FALSE),tip=rep("black", length(tree$tip.label)))         
    }
    ind = Descendants(tree,node,"all")
    x$edge$color[ind] = edge.color
    x$edge$width[ind] = edge.width
    x$edge$lty[ind] = edge.lty
    x[[2]][Descendants(tree, node, "tips")[[1]]] = tip.color
    if(plot){
        tree=reorder(tree)
        plot(tree, edge.color=x$edge$color[tree$edge[,2]], edge.width = x$edge$width[tree$edge[,2]], edge.lty = x$edge$lty[tree$edge[,2]], tip.color=x[[2]],...)
    }
    else return(x)
} 



