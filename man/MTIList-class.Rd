\name{MTIList-class}
\docType{class}
\alias{MTIList-class}
\alias{as.data.frame,MTIList-method}
\alias{entrezid,MTIList-method}
\alias{experiments,MTIList-method}
\alias{gene,MTIList-method}
\alias{geneSpecies,MTIList-method}
\alias{id,MTIList-method}
\alias{matmirna,MTIList-method}
\alias{mirfam,MTIList-method}
\alias{mirnaSpecies,MTIList-method}
\alias{pmid,MTIList-method}
\alias{premirna,MTIList-method}
\alias{reportCount,MTIList-method}
\alias{show,MTIList-method}
\alias{supportedBy,MTIList-method}
\alias{MTIList}

\title{Class \dQuote{MTIList}: list of miRNA target gene interaction}

\description{
  The class \dQuote{MTIList} extends the \code{SimpleList} class from
  the \code{S4Vectors} package and thus allows to store a collection of
  \code{\linkS4class{MTI}} objects. In addition, all functions
  applicable to a \code{list} can be used to this type of object.
}
\section{Objects from the Class}{
  Objects should be created using the function \code{MTIList}.
}
\section{Methods}{
  \describe{

    \item{\code{as.data.frame( x, collapse.reports=NULL,
	stringsAsFactors=getOption( "stringsAsFactors", TRUE ), ... )}}{
      Casts the \code{MTIList} object into a \code{data.frame}. If
      \code{collapse.reports} is not \code{NULL}, the values from the
      internal \code{Report} classes will be collapsed (separated by
      the value of \code{collapse.reports}) resulting in a
      \code{data.frame} with a single row for the \code{MTI}
      object. Otherwise, the \code{data.frame} will have as many rows
      as there are \code{Report} objects in the \code{MTI} object.
    }

    \item{\code{entrezid( object, ... )}}{
      Returns a character vector of the NCBI Entrezgene identifier of
      the target gene of each \code{\linkS4class{MTI}} object in the list.
    }

    \item{\code{experiments( object, ... )}}{
      Returns a list with the experiments that were performed to proof
      the interaction for each \code{\linkS4class{MTI}} object in the
      list.
    }

    \item{\code{gene( object, ... )}}{
      Returns a character vector with the official symbol of the target gene of each
      \code{\linkS4class{MTI}} object in the list.
    }

    \item{\code{geneSpecies( object, ... )}}{
      Returns a character vector with the species of the target gene of each
      \code{\linkS4class{MTI}} object in the list.
    }

    \item{\code{id( object, ... )}}{
      Returns a character vector with the identifier of each \code{\linkS4class{MTI}} object in
      the list.
    }

    \item{\code{matmirna( object, ... )}}{
      Returns a character vector with the name of the mature miRNA of each
      \code{\linkS4class{MTI}} object in the list.
    }

    \item{\code{mirnaSpecies( object, ... )}}{
      Returns a character vector with the species of the miRNA of each \code{\linkS4class{MTI}}
      object in the list.
    }

    \item{\code{pmid( object, ... )}}{
      Returns a list of the PubMed identifier(s) of the
      publications providing evidence of the interaction for each \code{\linkS4class{MTI}}
      object in the list.
    }

    \item{\code{reportCount( object, ... )}}{
      Returns a numeric vector with the number of publications/reports providing evidence for
      each interaction in the list.
    }

    \item{\code{show( object )}}{
      Print the content of the MTIList class.
    }

    \item{\code{supportedBy( object, ... )}}{
      Returns a list of the support type of the evidence/experiments for
      each \code{\linkS4class{MTI}} object in the list.
    }

    Retrieve annotations (from the \code{mirbase.db} package) for the mature
    miRNA.

    \item{\code{premirna( object, ... )}}{
      Returns a list with the the pre-miRNA name (i.e. the miRNA id) for the mature
      miRNA from the mirbase for each \code{\linkS4class{MTI}}
      object in the list.
    }

    \item{\code{mirfam( object, ... )}}{
      Returns a character vector of the miRNA family name for the mature miRNA from the
      mirbase for each \code{\linkS4class{MTI}}
      object in the list.
    }

  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Johannes Rainer}
\note{
  This object directly extends the \code{SimpleList} class from the
  \code{S4Vectors} package, thus any method applicable to that type of
  class (or the generic \code{list}) can be applied to objects of this type.
}
\seealso{
	\code{\linkS4class{MTI}}
}
\examples{

## get all MTIs from the database describing an interaction
## between the genes BCL2 and BCL2L11 and miRNAs from the mir-17 family.
BCL2L11 <- mtis( mirtarbase, filter=list( GenenameFilter( c( "BCL2", "BCL2L11" ) ), MirfamFilter( "mir-17" ) ) )
BCL2L11

## what's the mature miRNAs?
matmirna( BCL2L11 )
## what's are the precursors of these miRNA
premirna( BCL2L11 )

## the gene name and the Entrezgene ID?
gene( BCL2L11 )
entrezid( BCL2L11 )

## on how many reports (publications) does this base?
reportCount( BCL2L11 )

## what evidence/support types?
supportedBy( BCL2L11 )

## what experiments:
experiments( BCL2L11 )

## cast the MTIList object into a data.frame
as.data.frame( BCL2L11 )

## into a data.frame with a single row for each MTI
as.data.frame( BCL2L11, collapse.reports=";" )

}
\keyword{classes}

