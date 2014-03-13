\docType{methods}
\name{EstimateInsertSizeFromSame}
\alias{ContigLinkCount}
\alias{EstimateAssQual}
\alias{EstimateInsertSizeFromDiff}
\alias{EstimateInsertSizeFromSame}
\alias{plot}
\title{Estimate Insert size library from pairs of sequence mapping on the same contig.}
\usage{
\S4method{EstimateInsertSizeFromSame}{ContigLinkViewer}(object, ...)

\S4method{EstimateInsertSizeFromDiff}{ContigLinkViewer}(object, ...)

\S4method{EstimateAssQual}{ContigLinkViewer}(object, EstimateCutOff = TRUE,
  Min = NA, Max = NA, SdNumber = NA)

\S4method{ContigLinkCount}{ContigLinkViewer}(object, MinCount)

\S4method{plot}{ContigLinkViewer}(object, x = NA, ...)
}
\arguments{
  \item{object}{An object of class
  \code{\link{ContigLinkViewer}}}

  \item{An}{object of class \code{\link{ContigLinkViewer}}}

  \item{object}{An object of class
  \code{\link{ContigLinkViewer}}}

  \item{EstimateCutOff}{If TRUE, the interval for the
  insert size is [MeanIsize-(SdNumber*SdIsize),
  MeanIsize+(SdNumber*SdIsize)]}

  \item{Min}{If EstimateCutOff is FALSE The min value for
  the insert Size}

  \item{Max}{If EstimateCutOff is FALSE The max value for
  the insert Size}

  \item{SdNumber}{The number of sd to keep for computing
  the insert size interval.}

  \item{An}{object of class \code{\link{ContigLinkViewer}}}

  \item{MinCount}{The minimum number of pairs connecting
  two contigs}

  \item{object}{An object of class
  \code{\link{ContigLinkViewer}}}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }

  \item{x}{the coordinates of points in the plot.
  Alternatively, a single plotting structure, function or
  \emph{any \R object with a \code{plot} method} can be
  provided.}
}
\description{
Estimate Insert size library from pairs of sequence mapping
on the same contig.

Estimate Insert size library from pairs of sequence mapping
on two different contigs.

Estimate the assembly quality

Count the number of pairs of sequence that connect two
contigs.

Plot method for the class \code{\link{ContigLinkViewer}}
}
\details{
Estimate Insert size library from pairs of sequence mapping
on the same contig and having the good orientation. The
well mapped paired-end sequence (-> <-) belong to this two
categorism: \itemize{ \item{1} The strand of the first
sequence in pair is + and its sam flag is 99. The strand of
the second sequence in the pair is - and its flag is 147.
The insert size is obtained by substracting the start
position of the first sequence on the contig to the end
position of the second sequence on the contig and by adding
1. \item{2} The strand of the first sequence in pair is -
and its sam flag is 83. The strand of the second sequence
in the pair is + and its flag is 163. The insert size is
obtained by substracting the start position of the second
sequence on the contig to the end position of the first
sequence on the contig and by adding 1. }

Estimate Insert size library from pairs of sequence mapping
on two contigs. In this case the orientation of the
sequence on the contig does not matter. To obtain the
insert size, two cases: \itemize{ \item{The sequence is the
first in the pair} The insert size is the length of the
first contig plus the end position of the second sequence
on the second contig minus the start position of the first
sequence on the first contig minus 1.  \item{The sequence
is the second in the pair} The insert size is the length of
the first contig plus the end position of the first
sequence on the first contig minus the start position of
the second sequence on the second contig minus 1. }

This method creates the AssQual slot which contain 4
objects: \itemize{ \item{parameters} A vector giving the
insert size parameters which will be used to estimate the
assembly quality. Min and Max or Mean and sd.
\item{dfSameQualAssByContig} A data.frame summarizing the
repartition of the pairs according to the insert size range
and the orientation of the pairs.  This data.frame has 6
columns: \itemize{ \item{seqnames.1} A vector giving the
name of the contig \item{strand.1} A vector giving the
mapping strand of the first sequence in pair
\item{strand.2} A vector giving the mapping strand of the
second sequence in pair \item{outof} A vector giving the
number of pair that are outside the confidence interval for
the insert size \item{nonav} A vector giving the number of
pair for which the insert size can not be computed
\item{within} A vector giving the number of pair with an
insert size inside the confifence interval }
\item{dfSameQualAssTot} A data.frame which sum the
dfSameQualAssByContig by orientation only disregarding
contig information.  \item{dfDiffQualAssByContig} A
data.frame summarizing the repartition of the pairs
according to the insert size range and the orientation of
the contigs that are linked.  This data.frame has 6
columns: for the insert size \item{seqnames.1} A vector
giving the name of the first contig \item{seqnames.2} A
vector giving the name of the second contig \item{Ori1} A
vector giving the orientation of the first contig
\item{Ori2} A vector giving the orientation of the second
contig \item{outof} A vector giving the number of pair that
are outside the confidence interval for the insert size
\item{within} A vector giving the number of pair with an
insert size inside the confifence interval } To get the
quality, a confident interval has to be computed for the
insert size.

This method creates the CLC slot which is a list of 2
objects: \itemize{ \item{MinCount} A vector giving the
minimum number of pairs binding two contigs \item{CLC} A
data.frame with 5 columns giving the number of pairs above
the MinCount parameter that support the connection beween
two contigs.  \itemize{ \item{seqnames.1} A vector giving
the name of the first contig \item{seqnames.2} A vector
giving the name of the second contig \item{Ori1} A vector
giving the orientation of the first contig \item{Ori2} A
vector giving the orientation of the second contig
\item{within} A vector giving the number of pairs above the
MinCount parameter that support the connection beween the
two contigs.  } }

This plot function allow to generate a graph to vizualize
the contigs and their link using the Rgraphviz package.
Contigs are the nodes and the number of pairs supporting
the link between two are specified by the edge. The
ContigLinkViewer object has to be complete.
}
\author{
Delphine Charif & Joseph Tran

Delphine Charif & Joseph Tran

Delphine Charif & Joseph Tran

Delphine Charif & Joseph Tran

Delphine Charif & Joseph Tran
}

