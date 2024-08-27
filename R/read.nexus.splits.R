## @aliases read.nexus.splits write.nexus.splits write.splits
## read.nexus.networx write.nexus.networx

#' Function to import and export splits and networks
#'
#' \code{read.nexus.splits}, \code{write.nexus.splits},
#' \code{read.nexus.networx}, \code{write.nexus.networx}
#' can be used to import and export splits and networks with nexus format
#' and allow to exchange these object with other software like SplitsTree.
#' \code{write.splits} returns a human readable output.
#'
#' @param file a file name.
#' @param obj An object of class splits.
#' @param weights Edge weights.
#' @param taxa logical. If TRUE a taxa block is added
#' @param append logical. If TRUE the nexus blocks will be added to a file.
#' @param splits logical. If TRUE the nexus blocks will be added to a file.
#' @param x An object of class splits.
#' @param zero.print character which should be printed for zeros.
#' @param one.print character which should be printed for ones.
#' @param print.labels logical. If TRUE labels are printed.
#' @param \dots Further arguments passed to or from other methods.
#' @param labels names of taxa.
#' @return \code{write.nexus.splits} and \code{write.nexus.networx} write out
#' the \code{splits} and \code{networx} object to read with
#' other software like SplitsTree.
#' \code{read.nexus.splits} and \code{read.nexus.networx} return an
#' \code{splits} and \code{networx} object.
#' @note \code{read.nexus.splits} reads in the splits block of a nexus file. It
#' assumes that different co-variables are tab delimited and the bipartition
#' are separated with white-space. Comments in square brackets are ignored.
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @seealso \code{\link[ape]{prop.part}}, \code{\link{lento}},
#' \code{\link{as.splits}}, \code{\link{as.networx}}
#' @keywords cluster
#' @examples
#'
#' (sp <- as.splits(rtree(5)))
#' write.nexus.splits(sp)
#' spl <- allCircularSplits(5)
#' plot(as.networx(spl))
#' write.splits(spl, print.labels = FALSE)
#'
#' @rdname read.nexus.splits
#' @export
read.nexus.splits <- function(file) {
  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  semico <- grep(";", X)
  X <- gsub("\\[(.*?)\\]", "", X) # get rid of comments
  i1 <- grep("TAXLABELS", X, ignore.case = TRUE)
  sp <- grep("SPLITS;", X, ignore.case = TRUE)
  i1 <- i1[i1 < sp]
  if(length(i1)>1)i1 <- i1[length(i1)]
  taxlab <- ifelse(length(i1) > 0, TRUE, FALSE)
  if (taxlab) {
    end <- semico[semico >= i1][1]
    x <- X[(i1):end] # assumes not a 'new line' after "TRANSLATE"
    x <- gsub("TAXLABELS", "", x, ignore.case = TRUE)
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    x <- gsub("['\"]", "", x)
  }
  spEnd <- grep("END;", X, ignore.case = TRUE)
  spEnd <- spEnd[spEnd > sp][1]
  dims <- grep("DIMENSION", X, ignore.case = TRUE)
  cyc <- grep("CYCLE", X, ignore.case = TRUE)
  fcyc <- FALSE
  matr <- grep("MATRIX", X, ignore.case = TRUE)
  format <- grep("FORMAT", X, ignore.case = TRUE)
  start <- matr[matr > sp][1] + 1
  end <- semico[semico > start][1] - 1
  format <- format[(format > sp) & (format < spEnd)]

  res <- vector("list", end - start + 1)
  weights <- numeric(end - start + 1)
  j <- 1

  flab <- fwei <- fcon <- fint <- FALSE

  if (length(format) > 0) {
    tmp <- X[format]
    tmp <- gsub("\\;", "", tmp)
    tmp <- gsub("\\s+", "", tmp)
    flab <- grepl("labels=left", tmp, ignore.case = TRUE)
    fwei <- grepl("weights=yes", tmp, ignore.case = TRUE)
    fcon <- grepl("confidences=yes", tmp, ignore.case = TRUE)
    fint <- grepl("intervals=yes", tmp, ignore.case = TRUE)
    ind <- cumsum(c(flab, fwei, fcon, fint))
    mformat <- sum(c(flab, fwei, fcon, fint))
  }

  if (fint) intervals <- numeric(end - start + 1)
  if (fcon) confidences <- numeric(end - start + 1)
  if (flab) labels <- vector("character", end - start + 1)

  for (i in start:end) {
    tmp <- X[i]
    tmp <- sub("\\s+", "", tmp)
    tmp <- strsplit(tmp, "\t")[[1]]
    if (flab) {
      labels[j] <- gsub("'", "", tmp[ind[1]]) |> as.numeric()
    }
    if (fwei) weights[j] <- as.numeric(tmp[ind[2]])
    if (fcon) confidences[j] <- as.numeric(tmp[ind[3]])
    if (fint) intervals[j] <- as.numeric(tmp[ind[4]])
    tmp <- tmp[length(tmp)]
    tmp <- gsub("\\,", "", tmp)
    res[[j]] <- as.integer(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
    j <- j + 1
  }
  if (length(cyc) > 0) {
    tmp <- X[cyc]
    tmp <- gsub("\\;", "", tmp)
    tmp <- gsub("CYCLE", "", tmp, ignore.case = TRUE)
    tmp <- sub("\\s+", "", tmp)
    cyc <- as.integer(na.omit(as.numeric(strsplit(tmp, " ")[[1]])))
    fcyc <- TRUE
  }
  attr(res, "labels") <- x
  if (fwei) attr(res, "weights") <- weights
  if (fint) attr(res, "intervals") <- intervals
  if (fcon) attr(res, "confidences") <- confidences
  if (flab) attr(res, "splitlabels") <- labels
  if (fcyc) attr(res, "cycle") <- cyc
  class(res) <- "splits"
  res
}



#' @rdname read.nexus.splits
#' @export
write.nexus.splits <- function(obj, file = "", weights = NULL, taxa = TRUE,
                               append = FALSE) {
  taxa.labels <- attr(obj, "labels")
  ntaxa <- length(taxa.labels)
  obj <- ONEwise(obj)
  ind <- which(lengths(obj) == ntaxa)
  if (length(ind)) obj <- obj[-ind]
  nsplits <- length(obj)
  if (is.null(weights)) weight <- attr(obj, "weights")
  if (is.null(weight)) fwei <- FALSE
  else fwei <- TRUE
  if (!append) {
    cat("#NEXUS\n\n", file = file)
    cat("[Splits block for Spectronet or SplitsTree]\n", file = file,
      append = TRUE)
    cat("[generated by phangorn", packageDescription("phangorn",
      fields = "Version"), "]\n\n", file = file, append = TRUE)
  }
  # TAXON BLOCK
  if (taxa) {
    cat(paste("BEGIN TAXA;\n\tDIMENSIONS ntax=", ntaxa, ";\n",
      sep = ""), file = file, append = TRUE)
    cat("\tTAXLABELS", paste(taxa.labels, sep = " "), ";\nEND;\n\n",
      file = file, append = TRUE)
  }
  # SPLITS BLOCK
  cat(paste("BEGIN SPLITS;\n\tDIMENSIONS ntax=", ntaxa, " nsplits=", nsplits,
    ";\n", sep = ""), file = file, append = TRUE)
  format <- "\tFORMAT labels=left"
  if (fwei) format <- paste(format, "weights=yes")
  else format <- paste(format, "weights=no")
  fcon <- fint <- flab <- FALSE
  if (!is.null(attr(obj, "confidences"))) {
    format <- paste(format, "confidences=yes")
    fcon <- TRUE
    conf <- attr(obj, "confidences")
    if (any(is.na(conf))) {
      conf[is.na(conf)] <- 0
      attr(obj, "confidences") <- conf
    }
    if (storage.mode(conf) == "character") {
      conf[conf == ""] <- "0"
      attr(obj, "confidences") <- conf
    }
  }
  else format <- paste(format, "confidences=no")
  if (!is.null(attr(obj, "intervals"))) {
    format <- paste(format, "intervals=yes")
    fint <- TRUE
  }
  else format <- paste(format, "intervals=no")
  if (!is.null(attr(obj, "splitlabels"))) flab <- TRUE
  format <- paste(format, ";\n",  sep = "")
  cat(format, file = file, append = TRUE)
  if (!is.null(attr(obj, "cycle"))) {
    cycle <- paste(attr(obj, "cycle"), collapse = " ")
    cat("\tCYCLE\t", cycle, ";\n", sep = "", file = file, append = TRUE)
  }
  cat("\tMATRIX\n", file = file, append = TRUE)

  for (i in 1:nsplits) {
    slab <- ifelse(flab, attr(obj, "splitlabels")[i], i)
    swei <- ifelse(fwei, paste(weight[i], "\t"), "")
    scon <- ifelse(fcon, paste(attr(obj, "confidences")[i], "\t"), "")
    sint <- ifelse(fint, paste(attr(obj, "intervals")[i], "\t"), "")
    cat("\t\t", slab, "\t", swei, scon, sint, paste(obj[[i]], collapse = " "),
      ",\n", file = file, append = TRUE, sep = "")
  }
  cat("\t;\nEND;\n", file = file, append = TRUE)
}


#' @rdname read.nexus.splits
#' @export
write.nexus.networx <- function(obj, file = "", taxa = TRUE, splits = TRUE,
                                append = FALSE) {
  if (!append) {
    cat("#NEXUS\n\n", file = file)
    cat("[Network block for Spectronet or SplitsTree]\n", file = file,
      append = TRUE)
    cat("[generated by phangorn", packageDescription("phangorn",
      fields = "Version"), "]\n\n", file = file, append = TRUE)
  }
  ntaxa <- length(obj$tip.label)
  # TAXON BLOCK
  if (taxa) {
    cat(paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", ntaxa, ";\n",
      sep = ""), file = file, append = TRUE)
    if (splits) taxalabel <- attr(obj$splits, "labels")
    else taxalabel <- obj$tip.label
    cat("\tTAXLABELS", paste(taxalabel, sep = " "), ";\nEND;\n\n",
      file = file, append = TRUE)
  }
  # SPLITS BLOCK
  spl <- obj$splits
  if (splits) {
    write.nexus.splits(spl, file = file, weights = NULL, append = TRUE,
      taxa = FALSE)
  }
  nvertices <- max(obj$edge)

  #    if(is.null(attr(obj, "coords")))
  if (is.null(obj$.plot$vertices)) vertices <- coords.equal.angle(obj)
  else vertices <- obj$.plot$vertices

  # y-axis differs between R and SplitsTree
  vertices[, 2] <- -vertices[, 2]

  if (is.null(obj$.plot)) edge.col <- obj$.plot$edge.color
  else edge.col <- NULL
  nedges <- nrow(obj$edge)
  # NETWORK BLOCK
  cat(paste("BEGIN NETWORK;\nDIMENSIONS ntax=", ntaxa,
    "\tnvertices=", nvertices, "\tnedges=", nedges, ";\n", sep = ""),
  file = file, append = TRUE)
  cat("DRAW to_scale;\n", file = file, append = TRUE)
  cat("TRANSLATE\n", file = file, append = TRUE)
  if (is.null(obj$translate)) {
    for (i in 1:ntaxa) {
      cat(i, " ", obj$tip.label[i], ",\n", sep = "", file = file, append = TRUE)
    }
  }
  else {
    translate <- obj$translate
    for (i in seq_along(translate$label)) {
      cat(translate$node[i], " ", translate$label[i], ",\n", sep = "",
        file = file, append = TRUE)
    }
  }
  cat(";\nVERTICES\n", file = file, append = TRUE)
  for (i in 1:nvertices) {
    cat(i, "\t", vertices[i, 1], "\t", vertices[i, 2], ",\n", sep = "",
      file = file, append = TRUE)
  }
  if (!is.null(obj$tip.label)) {
    cat(";\nVLABELS\n", file = file, append = TRUE)
    if (is.null(obj$translate)) {
      for (i in 1:ntaxa) {
        cat(i, "\t", obj$tip.label[i], ",\n", sep = "", file = file,
          append = TRUE)
      }
    }
    else {
      for (i in seq_along(translate$node)) {
        cat(translate$node[i], " ", translate$label[i], ",\n", sep = "",
          file = file, append = TRUE)
      }
    }
  }
  # cnet$splitIndex if splits = TRUE
  cat(";\nEDGES\n", file = file, append = TRUE)

  if (is.null(obj$.plot$edge.color)) edge.col <- "black"
  else edge.col <- obj$.plot$edge.color
  if (length(edge.col) < nedges) edge.col <- rep(edge.col, length = nedges)

  splI <- TRUE
  if (is.null(obj$splitIndex)) splI <- FALSE
  for (i in 1:nedges) {
    ecoli <- edge.col[i]
    spInd <- ifelse(splI, paste("\ts=", obj$splitIndex[i], sep = ""), "")
    edgeCol <- ifelse(ecoli == "black", "", paste("\tfg=",
      paste(col2rgb(ecoli), collapse = " "), sep = ""))
    cat(i, "\t", obj$edge[i, 1], "\t", obj$edge[i, 2], spInd, edgeCol, ",\n",
      sep = "", file = file, append = TRUE)
  }
  cat(";\n", file = file, append = TRUE)
  cat("END;\n", file = file, append = TRUE)
  # force SplitsTree to accept the file
  cat("\nBEGIN st_Assumptions;\n    uptodate;\nEND; [st_Assumptions]\n",
    file = file, append = TRUE)
}


#' @rdname read.nexus.splits
#' @export
read.nexus.networx <- function(file, splits = TRUE) {
  spl <- NULL
  if (splits) spl <- read.nexus.splits(file)

  X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  semico <- grep(";", X)
  X <- gsub("\\[(.*?)\\]", "", X) # get rid of comments
  netStart <- grep("BEGIN NETWORK;", X, ignore.case = TRUE)
  if(length(netStart)==0){
    if(splits) {
      warning("File does not contain network block, return only splits!")
      return(spl)
    }
    else stop("File does not contain network block!")
  }
  netEnd <- grep("END;", X, ignore.case = TRUE)
  netEnd <- netEnd[netEnd > netStart][1]
  dims <- grep("DIMENSION", X, ignore.case = TRUE)
  dims <- dims[(dims > netStart) & (dims < netEnd)]

  ntaxa <- 0
  nvertices <- 0
  nedges <- 0

  if (length(dims) > 0) {
    tmp <- X[dims]
    tmp <- gsub("\\s+", "", tmp)

    ntaxa <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
      "\\3", tmp, perl = TRUE, ignore.case = TRUE))
    nvertices  <- as.numeric(sub("(.+?)(nvertices\\s*\\=\\s*)(\\d+)(.+)",
      "\\3", tmp, perl = TRUE, ignore.case = TRUE))
    nedges <- as.numeric(sub("(.+?)(nedges\\s*\\=\\s*)(\\d+)(.+)",
      "\\3", tmp, perl = TRUE, ignore.case = TRUE))
  }
  transl <- grep("TRANSLATE", X, ignore.case = TRUE)
  translation <- if (length(transl) == 1 && transl > netStart) TRUE
  else FALSE
  translate.nodes <- FALSE
  if (translation) {
    end <- semico[semico > transl][1]
    x <- X[(transl + 1):end]
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    x <- gsub("['\']", "", x)
    if (length(x) == 2 * ntaxa) {
      TRANS <- matrix(x, ncol = 2, byrow = TRUE)
      TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
      TRANS <- data.frame(node = as.integer(TRANS[, 1]), label = TRANS[, 2])
    }
    else {
      y <- as.integer(x)
      node <- integer(ntaxa)
      label <- character(ntaxa)
      k <- 1
      for (i in seq_along(x)) {
        if (!is.na(y[i])) tmp <- y[i]
        else {
          node[k] <- tmp
          label[k] <- x[i]
          k <- k + 1
        }
      }
      TRANS <- data.frame(node = node, label = label)
    }
  }
  vert <- grep("VERTICES", X, ignore.case = TRUE)
  start <- vert[vert > max(dims, netStart)][1] + 1
  end <- semico[semico > start][1] - 1
  VERT <- matrix(0, nvertices, 3, dimnames = list(NULL, c("id", "x", "y")))
  j <- 1
  for (i in start:end) {
    tmp <- X[i]
    tmp <- gsub("\\,", "", tmp)
    tmp <- strsplit(tmp, "[[:space:]]")[[1]]
    VERT[j, 1] <- as.numeric(tmp[1])
    VERT[j, 2] <- as.numeric(tmp[2])
    VERT[j, 3] <- as.numeric(tmp[3])
    j <- j + 1
  }

  edges <- grep("EDGES", X, ignore.case = TRUE)
  start <- edges[edges > max(dims, netStart)][1] + 1
  end <- semico[semico > start][1] - 1
  EDGE <- NULL
  if (splits) EDGE <- matrix(0L, nedges, 4, dimnames = list(NULL, c("id",
      "vert_id_2", "vert_id_2", "splits_id")))
  else EDGE <- matrix(0L, nedges, 3, dimnames = list(NULL, c("id", "vert_id_2",
      "vert_id_2")))
  j <- 1
  for (i in start:end) {
    tmp <- X[i]
    tmp <- gsub("\\,", "", tmp)
    tmp <- strsplit(tmp, "[[:space:]]")[[1]]
    EDGE[j, 1] <- as.integer(tmp[1])
    EDGE[j, 2] <- as.integer(tmp[2])
    EDGE[j, 3] <- as.integer(tmp[3])
    if (splits) {
      EDGE[j, 4] <- as.integer(sub("s=", "", tmp[4], ignore.case = TRUE))
    }
    j <- j + 1
  }

  swapEdge <- function(x, old, new) {
    x[x == new] <- -1L
    x[x == old] <- new
    x[x == -1L] <- old
    x
  }
  swapRow <- function(x, old, new) {
    tmp <- x[old, ]
    x[old, ] <- x[new, ]
    x[new, ] <- tmp
    x
  }
  splitIndex <- if (ncol(EDGE) == 4) EDGE[, 4]
  else NULL
  # quick and dirty
  el <- sqrt(rowSums( (VERT[EDGE[, 2], c(2:3)] - VERT[EDGE[, 3], c(2:3)])^2))
  edge <- EDGE[, c(2:3)]
  vert <- VERT[, c(2:3)]

  if(length(TRANS$label) > sum(tabulate(edge)==1L)){
#    ntip <- length(TRANS$label) - sum(tabulate(edge)==1L)
    if(!is.null(spl)){
      nspl <- length(spl)
      spl <- addTrivialSplits(spl)
      ind <- unlist(spl[(nspl+1):length(spl)])
      lab <- attr(spl, "labels")[ind]
      pos <- match(lab, TRANS$label)
    }
    else{stop("Problem")}
    new_edges <- max(edge) + seq_along(ind)
    edge <- rbind(edge, cbind(TRANS$node[pos], new_edges))
    vert <- rbind( vert, vert[TRANS$node[pos], ])
    if(!is.null(splitIndex))splitIndex <- c(splitIndex, (nspl+1):length(spl))
    TRANS$node[pos] <- new_edges
    el <- c(el, rep(0, length(new_edges)))
  }

  if (translate.nodes) {
    oldLabel <- as.integer(TRANS$node)
    for (i in 1:ntaxa) {
      edge <- swapEdge(edge, oldLabel[i], i)
      vert <- swapRow(vert, oldLabel[i], i)
    }
  }
  dimnames(edge) <- NULL
  # y-axis differs between in R and SplitsTree
  vert[, 2] <- -vert[, 2]
  obj <- list(edge = edge, tip.label = TRANS$label, edge.length = el,
    Nnode = max(edge) - ntaxa, splitIndex = splitIndex, splits = spl,
    translate = TRANS)
  obj$.plot <- list(vertices = vert, edge.color = "black", edge.width = 3,
    edge.lty = 1)
  class(obj) <- c("networx", "phylo")
  reorder(obj)
  obj
}


#' @rdname read.nexus.splits
#' @export
write.splits <- function(x, file = "", zero.print = ".", one.print = "|",
                         print.labels = TRUE, ...) {
  labels <- attr(x, "labels")
  cx <- as.matrix(x, zero.print = zero.print, one.print = one.print)
  w <- FALSE
  if (!is.null(attr(x, "names"))) {
    nam <- TRUE
    vnames <- format(attr(x, "names"))
  }
  nam <- FALSE
  if (!is.null(attr(x, "weights"))) {
    w <- TRUE
    weight <- format(attr(x, "weights"))
  }
  d <- FALSE
  if (!is.null(attr(x, "data"))) {
    d <- TRUE
    data <- attr(x, "data")
  }
  if (print.labels) {
    for (i in seq_along(labels)) cat(labels[i], "\n", file = file,
                                     append = TRUE)
  }
  if (w)
    cat("weight", "\t", file = file, append = TRUE)
  if (d)
    cat(paste(colnames(data), "\t"), file = file, append = TRUE)
  cat("\n", file = file, append = TRUE) # "Matrix",
  for (i in seq_along(x)) {
    if (nam)
      cat(vnames[i], "\t", file = file, append = TRUE)
    if (d)
      cat(paste(data[i, ], "\t"), file = file, append = TRUE)
    if (w)
      cat(weight[i], "\t", file = file)
    cat("\n", paste(cx[i, ], collapse = ""), "\n", file = file, append = TRUE)
  }
}
