\name{read.aa}
\alias{read.aa}
\title{Read Amino Acid Sequences in a File}
\usage{
read.aa(file, format = "interleaved", skip = 0,
         nlines = 0, comment.char = "#", seq.names = NULL)
}
\arguments{
  \item{file}{a file name specified by either a variable of mode character,
    or a double-quoted string.}
  \item{format}{a character string specifying the format of the DNA
    sequences. Three choices are possible: \code{"interleaved"},
    \code{"sequential"}, or \code{"fasta"}, or any unambiguous
    abbreviation of these.}
  \item{skip}{the number of lines of the input file to skip before
    beginning to read data.}
  \item{nlines}{the number of lines to be read (by default the file is
    read untill its end).}
  \item{comment.char}{a single character, the remaining of the line
    after this character is ignored.}
  \item{seq.names}{the names to give to each sequence; by default the
    names read in the file are used.}
}
\description{
  This function reads amino acid sequences in a file, and returns a matrix 
  list of DNA sequences with the names of the taxa read in the file as
  row names.}
\value{
  a matrix of amino acid sequences.
}
\references{
  Anonymous. FASTA format description.
  \url{http://www.ncbi.nlm.nih.gov/BLAST/fasta.html}

  Felsenstein, J. (1993) Phylip (Phylogeny Inference Package) version
  3.5c. Department of Genetics, University of Washington.
  \url{http://evolution.genetics.washington.edu/phylip/phylip.html}
}
\seealso{
  \code{\link[ape]{read.dna}}, \code{\link[ape]{read.GenBank}}, 
  \code{\link[phangorn]{phyDat}}, \code{\link[seqinr]{read.alignment}}
}
\author{Klaus Schliep \email{klaus.schliep@gmail.com}}
\keyword{IO}
