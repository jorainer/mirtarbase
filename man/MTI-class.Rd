\name{MTI-class}
\docType{class}
\alias{MTI-class}
\alias{as.data.frame,MTI-method}
\alias{entrezid,MTI-method}
\alias{experiments,MTI-method}
\alias{gene,MTI-method}
\alias{geneSpecies,MTI-method}
\alias{id,MTI-method}
\alias{initialize,MTI-method}
\alias{matmirna,MTI-method}
\alias{matmirnaId,MTI-method}
\alias{mirfam,MTI-method}
\alias{mirnaSpecies,MTI-method}
\alias{pmid,MTI-method}
\alias{premirna,MTI-method}
\alias{premirnaId,MTI-method}
\alias{reportCount,MTI-method}
\alias{reports,MTI-method}
\alias{reports<-,MTI,list-method}
\alias{reports<-,MTI,MTIReport-method}
\alias{show,MTI-method}
\alias{supportedBy,MTI-method}
\alias{entrezid}
\alias{gene}
\alias{geneSpecies}
\alias{id}
\alias{initialize}
\alias{matmirna}
\alias{matmirnaId}
\alias{mirfam}
\alias{mirnaSpecies}
\alias{premirna}
\alias{premirnaId}
\alias{reportCount}
\alias{reports}
\alias{reports<-}

\title{Class \dQuote{MTI}: miRNA target gene interaction}

\description{
  The class \dQuote{MTI} describes a miRNA target gene interaction as
  provided by the miRTarBase.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MTI", ...)} or
  using the \code{ \link{MTI} } function (the preferred way).
}
\section{Slots}{
  \describe{

    \item{\code{id}}{
      The miRTarBase identifier of the interaction.
    }

    \item{\code{mature_mirna}}{
      The (mirbase conform) mature miRNA
      identifier (e.g. \dQuote{hsa-mir-24-3p}).
    }

    \item{\code{query}}{
      Not used yet.
    }

    \item{\code{report}}{
      A list of \code{\linkS4class{MTIReport}} objects
      providing the evidence for the miRNA target gene interaction.
    }

    \item{\code{species_mirna}}{
      The species of the miRNA.
    }

    \item{\code{species_target_gene}}{
      The species of the target gene.
    }

    \item{\code{target_gene}}{
      The official gene symbol of the target gene.
    }

    \item{\code{target_gene_entrezid}}{
      The NCBI Entrezgene identifier of
      the target gene.
    }

  }
}
\section{Methods}{
  \describe{

    \item{\code{as.data.frame(x, collapse.reports=NULL,
	stringsAsFactors=getOption("stringsAsFactors", TRUE), ...)}}{
      Casts the \code{MTI} object into a \code{data.frame}. If
	\code{collapse.reports} is not \code{NULL}, the values from the
	internal \code{MTIReport} classes will be collapsed (separated by
	the value of \code{collapse.reports}) resulting in a
	\code{data.frame} with a single row for the \code{MTI}
	object. Otherwise, the \code{data.frame} will have as many rows
	as there are \code{MTIReport} objects in the \code{MTI} object.
    }

    \item{\code{entrezid(object, ...)}}{
      Returns the NCBI Entrezgene identifier of the target gene.
    }

    \item{\code{experiments(object, ...)}}{
      Returns the experiments
      that were performed to proof the interaction.
    }

    \item{\code{gene(object, ...)}}{
      Returns the official symbol of the target gene.
    }

    \item{\code{geneSpecies(object, ...)}}{
      Returns the species of the target gene.
    }

    \item{\code{id(object, ...)}}{
      Returns the identifier.
    }

    \item{\code{matmirna(object, ...)}}{
      Returns the name of the mature miRNA.
    }

    \item{\code{mirnaSpecies(object, ...)}}{
      Returns the species of the miRNA.
    }

    \item{\code{pmid(object, ...)}}{
      Returns the PubMed identifier(s) of the
      publications providing evidence of the interaction.
    }

    \item{\code{reportCount(object, ...)}}{
      Returns the number of publications/reports providing evidence for
      the interaction.
    }

    \item{\code{reports(x, ...)}, \code{reports(x) <- value}}{
      Getter and setter method of the \code{report} slot. Returns a list
      of \code{\linkS4class{MTIReport}} objects; \code{value} should be a
      list of \code{\linkS4class{MTIReport}} objects or a single object.
    }

    \item{\code{show(object)}}{
      Print the content of the MTI class.
    }

    \item{\code{supportedBy(object, ...)}}{
      Returns the support type of the evidence/experiments.
    }

    Retrieve annotations (from the \code{mirbase.db} package) for the mature
    miRNA.

    \item{\code{matmirnaId(object, ...)}}{
      Returns the identifier (accession) of the mature miRNA from the
      mirbase.
    }

    \item{\code{premirna(object, ...)}}{
      Returns the pre-miRNA name (i.e. the miRNA id) for the mature
      miRNA from the mirbase.
    }

    \item{\code{premirnaId(object, ...)}}{
      Returns the pre-miRNA accession for the mature miRNA from the
      mirbase.
    }

    \item{\code{mirfam(object, ...)}}{
      Returns the miRNA family name for the mature miRNA from the
      mirbase.
    }

  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Johannes Rainer}
\note{
  Instead of directly accessing the slots of the class, the getter methods listed above should be used.
}
\seealso{
	\code{\linkS4class{MTIReport}}
}
\examples{

## get an MTI from the database that represents an interaction
## between the gene BCL2L11 and miRNAs from the mir-17 family.
BCL2L11 <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2L11"), MirfamFilter("mir-17")))
BCL2L11 <- BCL2L11[[1]]
BCL2L11

## what's the mature miRNA?
matmirna(BCL2L11)
## what's the precursor of this miRNA
premirna(BCL2L11)

## the gene name and the Entrezgene ID?
gene(BCL2L11)
entrezid(BCL2L11)

## on how many reports (publications) does this base?
reportCount(BCL2L11)

## what evidence/support types?
supportedBy(BCL2L11)

## what experiments:
experiments(BCL2L11)

## list the reports:
reports(BCL2L11)

## cast the MTI object into a data.frame
as.data.frame(BCL2L11)

## into a single row data.frame
as.data.frame(BCL2L11, collapse.reports=";")

}
\keyword{classes}

