\docType{class}
\name{ContigLinkViewer-class}
\alias{ContigLinkViewer-class}
\alias{GetAssQual}
\alias{GetInsertSizeDiff}
\alias{GetInsertSizeSame}
\alias{SetAssQual}
\alias{SetContigLinkCount}
\alias{SetInsertSizeDiff}
\alias{SetInsertSizeSame}
\title{Class ContigLinkViewer}
\usage{
\S4method{GetAssQual}{ContigLinkViewer}(object)

\S4method{GetInsertSizeSame}{ContigLinkViewer}(object)

\S4method{GetInsertSizeDiff}{ContigLinkViewer}(object)

\S4method{initialize}{ContigLinkViewer}(.Object, BamFile,
  ExperimentalInsertSize)
}
\arguments{
  \item{BamFile}{A Bam file resulting from the mapping of
  paired-end sequences or revers complemented mated paire
  on a set of contigs (assembly experiment). The Bam file
  has to be sorted and indexed. Sequences that do not
  mapped or singletons have to be filtered out.}

  \item{ExperimentalInsertSize}{The experimental insert
  size. The mean length of dna fragment that have been
  keeped for the sequencing experiment.}

  \item{.Object}{ An object: see the Details section.}
}
\description{
Class for storing descriptive statistics on the pairs of
sequence that bind contigs as well as to estimate assembly
quality.

Setteur method for the slot ContigLinkCount of the class
\code{\link{ContigLinkViewer}}.

Setteur method for the slot InsertSizeSame of the class
\code{\link{ContigLinkViewer}}.

Setteur method for the slot InsertSizeDiff of the class
\code{\link{ContigLinkViewer}}.

Setteur method for the slot AssQual of the class
\code{\link{ContigLinkViewer}}.

Getteur method for the slot AssQual of the class
\code{\link{ContigLinkViewer}}.

Getteur method for the slot InsertSizeSame of the class
\code{\link{ContigLinkViewer}}.

Getteur method for the slot InsertSizediff of the class
\code{\link{ContigLinkViewer}}.

Constructor method for the class
\code{\link{ContigLinkViewer}}.
}
\details{
An assembly experiment from short pairs of sequences rarely
gives a uniq contig. This class aims (i) at giving some
confidence indices relying to an assembly solution by
computing some descriptive statistics (ii) at bringing out
links between contig to propose scaffold. In this class,
the information given from pairs mapping on the same contig
or from pairs binding two differents contig is treated
differently. Insert size estimated from sequencing pairs
mapping on the same contig is used as scaffolding
parameter.

This method instanciate the slots dfSame, dfDiff and
ExpInsertSize. For pairs mapping on different contigs,
stored in dfDiff, the Ori1 and Ori2 column are added. This
columns indicates the contigs orientation so has paired-end
reads are in good orientation (+ - or -> <-).
}
\author{
Delphine Charif & Joseph Tran

Delphine Charif & Joseph Tran
}

