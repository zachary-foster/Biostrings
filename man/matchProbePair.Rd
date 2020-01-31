\name{matchProbePair}

\alias{matchProbePair}
\alias{matchProbePair,DNAString-method}
\alias{matchProbePair,XStringViews-method}
\alias{matchProbePair,MaskedDNAString-method}


\title{Find "theoretical amplicons" mapped to a probe pair}

\description{
  In the context of a computer-simulated PCR experiment, one wants to find
  the amplicons mapped to a given primer pair.
  The \code{matchProbePair} function can be used for this: given a forward and a
  reverse probe (i.e. the chromosome-specific sequences of the forward and
  reverse primers used for the experiment) and a target sequence (generally a
  chromosome sequence), the \code{matchProbePair} function will return all
  the "theoretical amplicons" mapped to this probe pair.
}

\usage{
matchProbePair(Fprobe, Rprobe, subject, max.mismatch=0, min.mismatch=0, with.indels=FALSE,
               fixed=TRUE, algorithm="auto", logfile=NULL, verbose=FALSE)
}

\arguments{
  \item{Fprobe}{
    The forward probe.
  }
  \item{Rprobe}{
    The reverse probe.
  }
  \item{subject}{
    A \link{DNAString} object (or an \link{XStringViews} object
    with a \link{DNAString} subject) containing the target sequence.
  }
  \item{max.mismatch, min.mismatch}{
    The maximum and minimum number of mismatching letters allowed (see
    \code{?`\link{lowlevel-matching}`} for the details).
    If non-zero, an algorithm that supports inexact matching is used.
  }
  \item{with.indels}{
    If \code{TRUE} then indels are allowed. In that case, \code{min.mismatch}
    must be \code{0} and \code{max.mismatch} is interpreted as the maximum
    "edit distance" allowed between the pattern and a match.
    Note that in order to avoid pollution by redundant matches,
    only the "best local matches" are returned.
    Roughly speaking, a "best local match" is a match that is locally
    both the closest (to the pattern P) and the shortest.
    More precisely, a substring S' of the subject S is a "best local match" iff:
    \preformatted{
       (a) nedit(P, S') <= max.mismatch
       (b) for every substring S1 of S':
               nedit(P, S1) > nedit(P, S')
       (c) for every substring S2 of S that contains S':
               nedit(P, S2) >= nedit(P, S')
    }
    One nice property of "best local matches" is that their first and last
    letters are guaranteed to match the letters in P that they align with.
  }
  \item{fixed}{
    If \code{TRUE} (the default), an IUPAC ambiguity code in a probe sequence
    can only match the same code in the subject, and vice versa.
    If \code{FALSE}, an IUPAC ambiguity code in a probe sequence can match
    any letter in the subject that is associated with the code, and 
    vice versa. See \code{?`\link{lowlevel-matching}`} for more information.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"} or \code{"shift-or"}.
    See \code{\link{matchPattern}} for more information.
  }
  \item{logfile}{
    A file used for logging.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
}

\details{
  The \code{matchProbePair} function does the following: (1) find all
  the "plus hits" i.e. the Fprobe and Rprobe matches on the "plus" strand,
  (2) find all the "minus hits" i.e. the Fprobe and Rprobe matches on the
  "minus" strand and (3) from the set of all (plus\_hit, minus\_hit) pairs,
  extract and return the subset of "reduced matches" i.e. the (plus\_hit, minus\_hit)
  pairs such that (a) plus\_hit <= minus\_hit and (b) there are no hits (plus
  or minus) between plus\_hit and minus\_hit.
  This set of "reduced matches" is the set of "theoretical amplicons".
}

\value{
  An \link{XStringViews} object containing the set of "theoretical amplicons".
}

\author{H. Pagès}

\seealso{
  \code{\link{matchPattern}},
  \code{\link{matchLRPatterns}},
  \code{\link{findPalindromes}},
  \code{\link{reverseComplement}},
  \link{XStringViews-class}
}

\examples{
library(BSgenome.Dmelanogaster.UCSC.dm3)
subject <- Dmelanogaster$chr3R

## With 20-nucleotide forward and reverse probes:
Fprobe <- "AGCTCCGAGTTCCTGCAATA"
Rprobe <- "CGTTGTTCACAAATATGCGG"
matchProbePair(Fprobe, Rprobe, subject) # 1 "theoretical amplicon"

## With shorter forward and reverse probes, the risk of having multiple
## "theoretical amplicons" increases:
Fprobe <- "AGCTCCGAGTTCC"
Rprobe <- "CGTTGTTCACAA"
matchProbePair(Fprobe, Rprobe, subject) # 2 "theoretical amplicons"
Fprobe <- "AGCTCCGAGTT"
Rprobe <- "CGTTGTTCACA"
matchProbePair(Fprobe, Rprobe, subject) # 9 "theoretical amplicons"
}

\keyword{methods}

