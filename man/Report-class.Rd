\name{MTIReport-class}
\docType{class}
\alias{MTIReport-class}
\alias{show,MTIReport-method}
\alias{pmid,MTIReport-method}
\alias{as.data.frame,MTIReport-method}
\alias{experiments,MTIReport-method}
\alias{supportedBy,MTIReport-method}
\alias{pmid}
\alias{experiments}
\alias{supportedBy}

\title{Class \dQuote{MTIReport}, report and evidence for a miRNA target gene interaction }
\description{
  The class \dQuote{MTIReport} provides the evidence (i.e. the
  publication with all of the performed experiments of the interaction)
  for a miRNA target gene interaction as defined by the miRTarBase.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MTIReport", ...)} or
  using the \code{\link{MTIReport}} function (the preferred way).
}
\section{Slots}{
  \describe{
    \item{\code{pmid}}{
      The PubMed identifier of the publication in which
      the miRNA target gene interaction has been described.
    }
    \item{\code{experiments}}{
      A character vector of experiments that have been performed in the
      publication to verify the MTI.
    }
    \item{\code{support_type}}{
      The evidence grade of the interaction.
    }
  }
}
\section{Methods}{
  \describe{

    \item{\code{as.data.frame(x, stringsAsFactors=getOption(
	"stringsAsFactors", TRUE), ...)}}{
      Cast the MTIReport object into a (single row) \code{data.frame}.
    }

    \item{\code{show(object)}}{
      Print the content of the MTIReport class.
    }

    \item{\code{pmid(object, ...)}}{
      Returns the PubMed identifier(s) of the publications providing
      evidence of the interaction.
    }

    \item{\code{experiments(object, ...)}}{
      Returns the experiments that were performed to proof the
      interaction.
    }

    \item{\code{supportedBy(object, ...)}}{
      Returns the support type of the evidence/experiments.
    }
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Johannes Rainer}
\note{
  Instead of directly accessing the slots of the class, the getter methods listed above should be used.
}
\seealso{
	\code{\linkS4class{MTI}}
}
\examples{

## Get an MTI from the database we're looking here for an interaction
## between the gene BCL2L11 and miRNAs from the mir-17 family.
BCL2L11 <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2L11"), MirfamFilter("mir-17")))
BCL2L11 <- BCL2L11[[1]]
BCL2L11

## List the reports:
reports(BCL2L11)

## Get the first one
Rep <- reports(BCL2L11)[[1]]

Rep

## Get the support type:
supportedBy(Rep)

## Experiments
experiments(Rep)

## pubmed id
pmid(Rep)

## Cast the MTIReport object into a (single row) data.frame
as.data.frame(Rep)

}
\keyword{classes}

