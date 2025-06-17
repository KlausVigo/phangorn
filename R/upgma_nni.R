upgma_nni <- function(d, tree=NULL, opt = "ME", trace = 0, tipdated=FALSE,
         maxiter=Inf){
  me <- function(nh){
    res <- nh[c(5,5,6,6,4)] - nh[c(1,2,3,5,6)]
    if(any(res < 0)) return(+Inf)
    sum(res)
  }
  opt <- match.arg(opt, c("ME", "LS"))
  if(is.null(tree)) tree <- upgma(d)
  labels <- tree$tip.label
  nTips <- length(labels)
  if(inherits(d, "dist")){
    if(identical(attr(d, "Labels"), tree$tip.label)) y <- d
    else{
      y <- as.matrix(d)[labels, labels]
      y <- as.numeric(y[lower.tri(y)])
    }
  }
  best.tree <- tree
  if(opt=="LS") bestLS <- sum( (coph(best.tree) - y)^2)
  else bestME <- sum(best.tree$edge.length)
  run.nni <- TRUE
  iter <- 0L
  swap <- 0L
  if (trace > 0){
    if(opt=="ME") cat("start, ME:", bestME, "\n")
    else cat("start, RSS:", bestLS, "\n")
  }
  swap <- 0

  under_risk <- seq(Ntip(best.tree)+1L, Ntip(best.tree)+Nnode(best.tree))
  while (run.nni && iter < maxiter) {
    index <- indexNNI_fitch(best.tree, -1L)
    l <- nrow(index)
    INDEX <- rbind(index[, c(1, 3, 2, 4:6)], index[, c(2, 3, 1, 4:6)])
    nh <- nodeHeight(best.tree)
    if(!tipdated) nh[seq_len(nTips)] <- 0.0
    NH <- matrix(NA_real_, 2L*l, 6L)
    NH[] <- nh[INDEX]
    se0 <- numeric(l)
    se <- numeric(2*l)
    desc <- Descendants(best.tree)
    l_desc <- lengths(desc)
    # iterate only over active set
    active <- sort(unique(which(!is.na(match(index, under_risk))) %% nrow(index)))
    if(active[1]==0) active[1] <- nrow(index)
    for(i in active){
      se0[i] <- me(nh[index[i,]])
      #a, b, c
      l_ab <- l_desc[index[i,1]] * l_desc[index[i,2]]
      l_ac <- l_desc[index[i,1]] * l_desc[index[i,3]]
      l_bc <- l_desc[index[i,2]] * l_desc[index[i,3]]


      AB <- node_heights_upgma(desc[[index[i,1]]], desc[[index[i,2]]], nTips, y)
      AC <- node_heights_upgma(desc[[index[i,1]]], desc[[index[i,3]]], nTips, y)
      BC <- node_heights_upgma(desc[[index[i,2]]], desc[[index[i,3]]], nTips, y)
      # only check those under risk
      NH[i, 5] <- AC
      NH[i, 6] <-((AB * l_ab) + (BC * l_bc)) /
        ((l_desc[index[i,1]] + l_desc[index[i,3]]) * l_desc[index[i,2]])
      NH[i+l, 5] <- BC
      NH[i+l, 6] <-((AB * l_ab) + (AC * l_ac)) /
        ((l_desc[index[i,2]] + l_desc[index[i,3]]) * l_desc[index[i,1]])
      se[i] <- me(NH[i, ]) #
      if(opt=="LS" && is.finite(se[i])) {
        se[i] <- se[i]
      }
      se[i+l] <- me(NH[i+l, ])
      if(opt == "LS" && is.finite(se[i+l])) se[i+l] <- se[i+l]
    }
    M <- se - c(se0, se0)
    candidates <- which(M < -1e-8)
    if(length(candidates) == 0) run.nni <- FALSE
    # nur 5 & 6 sollten überschneiden
    under_risk <- INDEX[candidates, c(5,6)] |> as.vector() |> unique()
    while (length(candidates) > 0) {
      pscore <- M[candidates]
      ind <- candidates[which.min(pscore)]
      tree2 <- changeEdge(best.tree, INDEX[ind, c(2, 3)])
      nh[ INDEX[ind, c(5, 6)] ] <- NH[ind, c(5,6)]
      tree2$edge.length <- nh[tree2$edge[,1]] - nh[tree2$edge[,2]]
      swap <- swap + 1
      best.tree <- tree2
      indi <- c( which(INDEX == INDEX[ind,5], arr.ind = TRUE),
                 which(INDEX == INDEX[ind,6], arr.ind = TRUE)) |> unique()
      candidates <- setdiff(candidates, indi)
    }
    if (trace > 0) cat("NNI rearangements:", swap, "\n")
    iter <- iter + 1L
  }
  attr(best.tree, "swap" ) <- swap
  best.tree
}
#' @srrstats {G2.3, G2.3a} in lines: 8
