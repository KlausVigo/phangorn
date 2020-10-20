distC <- function(d, CL) {
  l <- length(CL)
  res <- matrix(0, l, l)
  for (i in 1:(l - 1)) {
    for (j in (i + 1):l)
      res[i, j] <- mean.default(d[CL[[i]], CL[[j]]])
  }
  res + t(res)
}


updateDM <- function(DM, d, CL, j) {
  l <- length(CL)
  for (i in 1:l) {
    DM[i, j] <- DM[j, i] <- mean.default(d[CL[[i]], CL[[j]]])
  }
  DM[j, j] <- 0
  DM
}


Rx <- function(d, x, CL) {
  lx <- length(x)
  res <- numeric(lx)
  lC <- length(CL)
  for (i in 1:lx) {
    xi <- x[i]
    tmp <- 0
    for (j in 1:lx) {
      if (j != i) tmp <- tmp + d[xi, x[j]]
    }
    if (lC > 0) {
      for (j in 1:lC) {
        tmp <- tmp + mean.default(d[xi, CL[[j]]])
      }
    }
    res[i] <- tmp
  }
  res
}


reduc <- function(d, x, y, z) {
  u <- 2 / 3 * d[x, ] + d[y, ] / 3
  v <- 2 / 3 * d[z, ] + d[y, ] / 3
  uv <- (d[x, y] + d[x, z] + d[y, z]) / 3
  d[x, ] <- u
  d[, x] <- u
  d[z, ] <- v
  d[, z] <- v

  d[y, ] <- 0
  d[, y] <- 0

  d[x, z] <- d[z, x] <- uv
  d[x, x] <- d[z, z] <- 0
  #  diag(d) <- 0
  d
}


# computes ordering
getOrderingNN <- function(x) {
  x <- as.matrix(x)
  labels <- attr(x, "Labels")
  if (is.null(labels))
    labels <- colnames(x)
  d <- x # as.matrix(x)
  l <- dim(d)[1]
  CL <- vector("list", l)
  CL[1:l] <- ORD <- 1:l
  lCL <- length(CL)
  ord <- CL

  DM <- DM_C <- DM_V <- d

  while (lCL > 1) {
#    i <- 0
#    j <- 0
    DM <- distC(d, CL)

    l <- nrow(DM)
    if (l > 2) {
      r <- rowSums(DM) / (l - 2)
      tmp <- out_cpp(DM, r, l)
#      tmp <- .C("out", as.double(DM), as.double(r), as.integer(l),
#        as.integer(i), as.integer(j), PACKAGE = "phangorn")
      e1 <- tmp[1]
      e2 <- tmp[2]
    }
    else {
      e1 <- 1
      e2 <- 2
    }
    n1 <- length(CL[[e1]])
    n2 <- length(CL[[e2]])
    if (n1 == 1 & n2 == 1) {
      newCL <- c(CL[[e1]], CL[[e2]])
      newOrd <- newCL
      CL <- c(CL[-c(e1, e2)], list(newCL))
      ord <- c(ord[-c(e1, e2)], list(newCL))
      lCL <- lCL - 1L
    }
    else {
      CLtmp <- c(as.list(CL[[e1]]), as.list(CL[[e2]]), CL[-c(e1, e2)])
      ltmp <- length(CLtmp)
      DM2 <- distC(d, CLtmp)
      if (ltmp > 2) rtmp <- rowSums(DM2) / (ltmp - 2)
      DM2 <- DM2 - outer(rtmp, rtmp, "+")
      # compute only this
      TMP <- DM2[1:n1, (n1 + 1):(n1 + n2)]

      blub <- which.min(TMP)
      # print("blub")
      if (n1 == 2 & n2 == 1) {
        if (blub == 2) {
          newCL <- c(CL[[e1]][1], CL[[e2]])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]])
        }
        else {
          newCL <- c(CL[[e2]], CL[[e1]][2])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]], CL[[e1]][1], CL[[e1]][2])
        }


      }
      if (n1 == 1 & n2 == 2) {
        if (blub == 1) {
          newCL <- c(CL[[e1]], CL[[e2]][2])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]], CL[[e2]][1], CL[[e2]][2])
        }
        else {
          newCL <- c(CL[[e2]][1], CL[[e1]])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]][1], CL[[e2]][2], CL[[e1]])
        }
      }
      if (n1 == 2 & n2 == 2) {
        if (blub == 1) {
          newCL <- c(CL[[e1]][2], CL[[e2]][2])
          newOrd <-  c(rev(ord[[e1]]), ord[[e2]])
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][1])
          d <- reduc(d, CL[[e1]][2], CL[[e2]][1], CL[[e2]][2])
        }
        if (blub == 2) {
          newCL <- c(CL[[e1]][1], CL[[e2]][2])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][1])
          d <- reduc(d, CL[[e1]][1], CL[[e2]][1], CL[[e2]][2])
        }
        if (blub == 3) {
          newCL <- c(CL[[e1]][2], CL[[e2]][1])
          newOrd <-  c(rev(ord[[e1]]), rev(ord[[e2]]))
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][2])
          d <- reduc(d, CL[[e1]][2], CL[[e2]][2], CL[[e2]][1])
        }
        if (blub == 4) {
          newCL <- c(CL[[e1]][1], CL[[e2]][1])
          newOrd <-  c(ord[[e1]], rev(ord[[e2]]))
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][2])
          d <- reduc(d, CL[[e1]][1], CL[[e2]][2], CL[[e2]][1])
        }
      }
      CL <- c(CL[-c(e1, e2)], list(newCL))
      ord <- c(ord[-c(e1, e2)], list(newOrd))
      lCL <- lCL - 1L
    }
  }
  newOrd
}


# computes ordering O(n^2) statt O(n^3) !!!
# needs debugging
getOrderingNN2 <- function(x) {
  x <- as.matrix(x)
  labels <- attr(x, "Labels")
  if (is.null(labels))
    labels <- colnames(x)
  d <- x # as.matrix(x)
  l <- dim(d)[1]
  CL <- vector("list", l)
  CL[] <- seq_len(l)
  lCL <- length(CL)
  ord <- CL
  # DM_C connected components, DM_V vertices
  DM_C <- DM_V <- DM <- d
  z <- 0
  while (lCL > 1) {
#    i <- 0
#    j <- 0
    z <- z + 1

    l <- nrow(DM)
    # compute Q_D from D_C
    if (l > 2) {
      r <- rowSums(DM) / (l - 2)
      tmp <- out_cpp(DM, r, l)
#      tmp <- .C("out", as.double(DM), as.double(r), as.integer(l),
#        as.integer(i), as.integer(j), PACKAGE = "phangorn")
      e1 <- tmp[1]
      e2 <- tmp[2]
    }
    else {
      e1 <- 1
      e2 <- 2
    }
    n1 <- length(CL[[e1]])
    n2 <- length(CL[[e2]])
    if (n1 == 1 & n2 == 1) {
      # add edge
      newCL <- c(CL[[e1]], CL[[e2]])
      newOrd <- newCL
      CL[[e1]] <- newCL
      # update DM_C
      DM <- updateDM(DM, d, CL, e1)
      DM <- DM[-e2, -e2, drop = FALSE]
      CL <- CL[-e2]

      ord[[e1]] <- newCL
      ord <- ord[-e2]

      lCL <- lCL - 1L
    }
    else {
      CLtmp <- c(as.list(CL[[e1]]), as.list(CL[[e2]]), CL[-c(e1, e2)])
      ltmp <- length(CLtmp)

      CLtmp2 <- c(CL[[e1]], CL[[e2]])
      rtmp2 <- Rx(d, CLtmp2, CL[-c(e1, e2)])
      if (ltmp > 2) rtmp2 <- rtmp2 / (ltmp - 2)
      DM3 <- d[CLtmp2, CLtmp2] - outer(rtmp2, rtmp2, "+")
      TMP2 <- DM3[1:n1, (n1 + 1):(n1 + n2)]
      blub <- which.min(TMP2)

      if (n1 == 2 & n2 == 1) {
        if (blub == 2) {
          newCL <- c(CL[[e1]][1], CL[[e2]])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]])
        }
        else {
          newCL <- c(CL[[e2]], CL[[e1]][2])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]], CL[[e1]][1], CL[[e1]][2])
        }

      }
      if (n1 == 1 & n2 == 2) {
        if (blub == 1) {
          newCL <- c(CL[[e1]], CL[[e2]][2])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]], CL[[e2]][1], CL[[e2]][2])
        }
        else {
          newCL <- c(CL[[e2]][1], CL[[e1]])
          newOrd <- c(ord[[e2]], ord[[e1]])
          d <- reduc(d, CL[[e2]][1], CL[[e2]][2], CL[[e1]])
        }
      }
      if (n1 == 2 & n2 == 2) {
        if (blub == 1) {
          newCL <- c(CL[[e1]][2], CL[[e2]][2])
          newOrd <-  c(rev(ord[[e1]]), ord[[e2]])
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][1])
          d <- reduc(d, CL[[e1]][2], CL[[e2]][1], CL[[e2]][2])
        }
        if (blub == 2) {
          newCL <- c(CL[[e1]][1], CL[[e2]][2])
          newOrd <-  c(ord[[e1]], ord[[e2]])
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][1])
          d <- reduc(d, CL[[e1]][1], CL[[e2]][1], CL[[e2]][2])
        }
        if (blub == 3) {
          newCL <- c(CL[[e1]][2], CL[[e2]][1])
          newOrd <-  c(rev(ord[[e1]]), rev(ord[[e2]]))
          d <- reduc(d, CL[[e1]][2], CL[[e1]][1], CL[[e2]][2])
          d <- reduc(d, CL[[e1]][2], CL[[e2]][2], CL[[e2]][1])
        }
        if (blub == 4) {
          newCL <- c(CL[[e1]][1], CL[[e2]][1])
          newOrd <-  c(ord[[e1]], rev(ord[[e2]]))
          d <- reduc(d, CL[[e1]][1], CL[[e1]][2], CL[[e2]][2])
          d <- reduc(d, CL[[e1]][1], CL[[e2]][2], CL[[e2]][1])
        }
      }
      ord[[e1]] <- newOrd
      ord <- ord[-e2]

      CL[[e1]] <- newCL

      DM <- updateDM(DM, d, CL, e1)
      DM <- DM[-e2, -e2, drop = FALSE]

      CL <- CL[-e2]
      lCL <- lCL - 1L
    }
  }
  newOrd
}


#' Computes a neighborNet from a distance matrix
#'
#' Computes a neighborNet, i.e. an object of class \code{networx} from a
#' distance matrix.
#'
#' \code{neighborNet} is still experimental. The cyclic ordering sometimes
#' differ from the SplitsTree implementation, the \emph{ord} argument can be
#' used to enforce a certain circular ordering.
#'
#' @param x a distance matrix.
#' @param ord a circular ordering.
#' @return \code{neighborNet} returns an object of class networx.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link{splitsNetwork}}, \code{\link{consensusNet}},
#' \code{\link{plot.networx}}, \code{\link{lento}},
#' \code{\link{cophenetic.networx}}, \code{\link{distanceHadamard}}
#' @references Bryant, D. & Moulton, V. (2004) Neighbor-Net: An Agglomerative
#' Method for the Construction of Phylogenetic Networks. \emph{Molecular
#' Biology and Evolution}, 2004, \bold{21}, 255-265
#' @keywords hplot
#' @examples
#'
#' data(yeast)
#' dm <- dist.ml(yeast)
#' nnet <- neighborNet(dm)
#' plot(nnet, "2D")
#'
#' @export neighborNet
neighborNet <-  function(x, ord = NULL) {
  x <- as.matrix(x)
  labels <- attr(x, "Labels")[[1]]
  if (is.null(labels))
    labels <- colnames(x)
  l <- length(labels)
  if (is.null(ord)) ord <- getOrderingNN(x)
  spl <- allCircularSplits(l, labels[ord])
  spl <- nnls.splits(spl, x)
  # nnls.split mit nnls statt quadprog
  attr(spl, "cycle") <- 1:l
  as.networx(spl)
}


removeNonsense <- function(obj) {
  nTips <- length(attr(obj, "label"))
  l <- lengths(obj)
  ind <- which( (l == 0L) | (l == nTips))
  obj <- obj[-ind]
}

