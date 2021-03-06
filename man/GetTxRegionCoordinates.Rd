% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{GetTxRegionCoordinates}
\alias{GetTxRegionCoordinates}
\title{Calculate transcript region starting coordinates.}
\usage{
GetTxRegionCoordinates(refGenome)
}
\arguments{
\item{refGenome}{A \code{character} string.}
}
\value{
A \code{list}. See 'Details'.
}
\description{
Calculate transcript region starting coordinates along full transcript from
a reference transcriptome. See 'Details'.
}
\details{
The function loads transcriptome data generated from \code{BuildTx()} and
stored in a \code{.RData} file, and calculates starting coordinates of the
transcript regions 5'UTR, CDS, 3'UTR along the full transcript. There should 
not be any need to call this function directly. The function returns a 
\code{list} with the following elements:
\enumerate{
    \item \code{start_region}: A \code{data.frame} of transcript region 
    start coordinates
    \item \code{gr.start}: A \code{GRanges} object of start codon transcript
    coordinates
    \item \code{gr.stop}: A \code{GRanges} object of stop codon transcript
    coordinates
}
}
\author{
Maurits Evers, \email{maurits.evers@anu.edu.au}
}
\keyword{internal}
