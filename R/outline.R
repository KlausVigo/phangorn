#spl <- structure(list(1L, 2L, 3L, 4L, 5L, 2:3,  3:5, 3:4),
#                 labels = c("t1", "t2", "t3", "t4", "t5"),
#                 weights = c(1, 1, 1, 1, 1, 1, 1, 1),
#                 class = c("splits", "prop.part"))
#  x <- spl

#' @importFrom stats complete.cases

outline <- function(x){
  l <- length(attr(x, "labels"))
  ord <- 1:l
  if (!is.null(attr(x, "cycle"))) ord <- attr(x, "cycle")
  x <- changeOrder(x, attr(x, "labels")[ord])
  x <- ONEwise(x)

  fun1 <- function(x){
    nTips <- length(attr(x, "labels"))
    v <- seq_len(nTips)
    x <- ONEwise(x)
    for (i in seq_along(x)) x[[i]] <- v[-x[[i]]]
    x
  }
  fun2 <- function(x){
    res <- integer(length(x))
    for (i in seq_along(x)) res[i] <- x[[i]][1]
    res
  }
  fun3 <- function(x){
    res <- integer(length(x))
    for (i in seq_along(x)) res[i] <- x[[i]][length(x[[i]])]
    res
  }

  y <- fun1(x)
  i <- fun2(y)
  j <- fun3(y)

  ord1 <- order(i, -j)
  ord2 <- order(j, -i)

  l <- length(y)
  n <- 1L
  m <- 1L
  k <- 1L
  outbound <- logical(2*l)
  final_ord <- integer(2*l)
  pos <- integer(2*l)
  for (k in seq_len(2*l)) {
    if ( (i[ord1[m]] < (j[ord2[n]] + 1))  && (m < l + 1) ) {
      final_ord[k] <- ord1[m]
      outbound[k] <- TRUE
      pos[k] <- m
      m <- m + 1L
    } else{
      final_ord[k] <- ord2[n]
      pos[k] <- n
      n <- n + 1L
    }
  }
  tmp <- data.frame(final_ord = final_ord, outbound = outbound)
## compute edge matrix and coord matrix

  nTips <- length(attr(y, "labels"))
#  angle <- (i + j - 2) / (2 * nTips) * (2 * pi)
  angle <- spl2angle(y)
  weight <- attr(y, "weight")
  if(is.null(weight)) weight <- rep(1, length(y))

  edge <- matrix(NA_integer_, 2*l-1, 2)
  edge_length <- rep(NA_integer_, 2*l-1)
  coord <- matrix(NA_real_,  2*l + 1, 2)

  nTips <- length(attr(y, "labels"))
  current_tip <- 1L
  current_node <- 1L
  max_node <- 1L

  coord[1, ] <- c(0,0)
  spl_list <- list(NULL)

  for (k in 1:(length(final_ord) - 1)) {
    p <- final_ord[k]
    if (outbound[k]) {
      edge[k,] <- c(current_node, max_node + 1L)
      coord[max_node + 1L, ] <- coord[current_node, ] +
        kreis2kart(weight[p], angle[p])
      spl_list[[max_node + 1]] <- sort(c(spl_list[[current_node]], p))
      max_node <- max_node + 1L
      current_node <- max_node
   } else {
      tmp <- spl_list[[current_node]]
      tmp <- tmp[tmp != p]
      ind <- match(list(tmp), spl_list)
      if(is.na(ind)){
        edge[k,] <- c(current_node, max_node + 1L)
        coord[max_node + 1L, ] <- coord[current_node, ] +
          kreis2kart(weight[p], angle[p] + pi)
        spl_list[[max_node + 1]] <- sort(tmp)
        max_node <- max_node + 1L
        current_node <- max_node
      } else{
        edge[k,] <- c(ind, current_node)
        current_node <- ind
      }
    }
  }
  coord <- coord[complete.cases(coord),]
  final_ord <- final_ord[seq_len(nrow(edge))]
  splitIndex <- final_ord[!duplicated(edge)]
  edge <- edge[complete.cases(edge),] |> unique()

  tab <- tabulate(edge, max(edge))
  pos <- integer(max(edge))
  pos[tab == 1L] <- seq_len(nTips)
  pos[tab > 1L] <- (nTips + 1L):max(edge)
  edge[] <- pos[edge]
  ind <- which(edge[,1] %in% seq_len(nTips))
  if(length(ind) > 0) edge[ind, ] <- edge[ind, c(2,1)]

  coord <- coord[order(pos),]

  res <- list(edge = edge, tip.label = attr(x, "labels"),
              Nnode = max(edge) - nTips,
              edge.length = NULL, splits = x, splitIndex = splitIndex)
  res$.plot <- list(vertices = coord)
  class(res) <- c("networx", "phylo")
  res
}
