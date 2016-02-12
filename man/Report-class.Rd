\name{Report-class}
\docType{class}
\alias{Report-class}
\alias{show,Report-method}
\alias{pmid,Report-method}
\alias{as.data.frame,Report-method}
\alias{experiments,Report-method}
\alias{supportedBy,Report-method}
\alias{pmid}
\alias{experiments}
\alias{supportedBy}

\title{Class \dQuote{Report}, report and evidence for a miRNA target gene interaction }
\description{
  The class \dQuote{Report} provides the evidence for a miRNA target gene interaction as
  defined by the miRTarBase.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Report", ...)} or
using the \code{ \link{Report} } function (the preferred way).
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

    \item{\code{as.data.frame( x, stringsAsFactors=getOption(
	"stringsAsFactors", TRUE ), ... )}}{
      Cast the Report object into a (single row) \code{data.frame}.
    }

    \item{\code{show( object )}}{
      Print the content of the Report class.
    }

    \item{\code{pmid( object, ... )}}{
      Returns the PubMed identifier(s) of the publications providing
      evidence of the interaction.
    }

    \item{\code{experiments( object, ... )}}{
      Returns the experiments that were performed to proof the
      interaction.
    }

    \item{\code{supportedBy( object, ... )}}{
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

## get an MTI from the database we're looking here for an interaction
## between the gene BCL2L11 and miRNAs from the mir-17 family.
BCL2L11 <- mtis( mirtarbase, filter=list( GenenameFilter( "BCL2L11" ), MirfamFilter( "mir-17" ) ) )
BCL2L11 <- BCL2L11[[1]]
BCL2L11

## list the reports:
reports( BCL2L11 )

## get the first one
Rep <- reports( BCL2L11 )[[1]]

Rep

## get the support type:
supportedBy( Rep )

## experiments
experiments( Rep )

## pubmed id
pmid( Rep )

## cast the Report object into a (single row) data.frame
as.data.frame( Rep )

}
\keyword{classes}

