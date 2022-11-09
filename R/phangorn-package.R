#' Chloroplast alignment
#'
#' Amino acid alignment of 15 genes of 19 different chloroplast.
#'
#'
#' @name chloroplast
#' @docType data
#' @keywords datasets
#' @examples
#'
#' data(chloroplast)
#' chloroplast
#'
NULL



#' Laurasiatherian data (AWCMEE)
#'
#' Laurasiatherian RNA sequence data
#'
#'
#' @name Laurasiatherian
#' @docType data
#' @source Data have been taken from the former repository of the Allan Wilson
#' Centre and were converted to R format by \email{klaus.schliep@gmail.com}.
#' @keywords datasets
#' @examples
#'
#' data(Laurasiatherian)
#' str(Laurasiatherian)
#'
NULL


#' @keywords internal
"_PACKAGE"


#' Internal phangorn Functions
#'
#' Internal \pkg{phangorn} functions.
#'
#' @name phangorn-internal
#' @aliases phangorn-internal threshStateC coords map_duplicates
#' @keywords internal
NULL

#' Yeast alignment (Rokas et al.)
#'
#' Alignment of 106 genes of 8 different species of yeast.
#'
#'
#' @name yeast
#' @docType data
#' @references Rokas, A., Williams, B. L., King, N., and Carroll, S. B. (2003)
#' Genome-scale approaches to resolving incongruence in molecular phylogenies.
#' \emph{Nature}, \bold{425}(6960): 798--804
#' @keywords datasets
#' @examples
#'
#' data(yeast)
#' str(yeast)
#'
NULL



#' Morphological characters of Mites (Schäffer et al. 2010)
#'
#' Matrix for morphological characters and character states for 12 species of
#' mites. See vignette '02_Phylogenetic trees from morphological data' for
#' examples to import morphological data.
#'
#' @name mites
#' @docType data
#' @references Schäffer, S., Pfingstl, T., Koblmüller, S., Winkler, K. A.,
#' Sturmbauer, C., & Krisper, G. (2010). Phylogenetic analysis of European
#' Scutovertex mites (Acari, Oribatida, Scutoverticidae) reveals paraphyly and
#' cryptic diversity: a molecular genetic and morphological approach.
#' \emph{Molecular Phylogenetics and Evolution}, \bold{55(2)}, 677--688.
#' @keywords datasets
#' @examples
#' data(mites)
#' mites
#' # infer all maximum parsimony trees
#' trees <- bab(mites)
#' # For larger data sets you might use pratchet instead bab
#' # trees <- pratchet(mites, minit=200, trace=0, all=TRUE)
#' # build consensus tree
#' ctree <- root(consensus(trees, p=.5), outgroup = "C._cymba",
#'               resolve.root=TRUE, edgelabel=TRUE)
#' plotBS(ctree, trees)
#' cnet <- consensusNet(trees)
#' plot(cnet)
NULL
