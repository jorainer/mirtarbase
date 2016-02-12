\name{MirtarbaseDb-class}
\Rdversion{1.1}
\docType{class}
\alias{MirtarbaseDb-class}
\alias{dbconn}
\alias{dbconn,MirtarbaseDb-method}
\alias{listColumns}
\alias{listColumns,MirtarbaseDb-method}
\alias{listExperiments}
\alias{listExperiments,MirtarbaseDb-method}
\alias{listPmids}
\alias{listPmids,MirtarbaseDb-method}
\alias{listSpecies}
\alias{listSpecies,MirtarbaseDb-method}
\alias{listSupportTypes}
\alias{listSupportTypes,MirtarbaseDb-method}
\alias{mtis}
\alias{mtis,MirtarbaseDb-method}
\alias{mtisBy}
\alias{mtisBy,MirtarbaseDb-method}
\alias{show}
\alias{show,MirtarbaseDb-method}
\alias{listTables}
\alias{listTables,MirtarbaseDb-method}
\alias{version}
\alias{version,MirtarbaseDb-method}

\title{miRNA target-gene interaction database}
\description{
  Use and retrieve miRNA target gene interactions (MTIs) as defined by
  the miRTarbase database.
}
\section{Objects from the Class}{
  The \code{MirtarbaseDb} object providing the dbconn to the
  respective MTI database is created upon loading of the package and is
  associated with the \code{mirtarbase} variable.
  See examples below for usage.
}
\usage{

\S4method{dbconn}{MirtarbaseDb}(x, ...)

\S4method{listColumns}{MirtarbaseDb}(x, table, skip.keys=TRUE, ...)

\S4method{listExperiments}{MirtarbaseDb}(x, ...)

\S4method{listPmids}{MirtarbaseDb}(x, ...)

\S4method{listSpecies}{MirtarbaseDb}(x, of="gene", ...)

\S4method{listSupportTypes}{MirtarbaseDb}(x, ...)

\S4method{mtis}{MirtarbaseDb}(x, attrs=listColumns(x, "mirtarbase"
), filter, order.by="mirtarbase_id", order.type="asc",
return.type="MTIList", force=FALSE)

\S4method{mtisBy}{MirtarbaseDb}(x, by="gene", filter)

\S4method{listTables}{MirtarbaseDb}(x, ...)

\S4method{version}{MirtarbaseDb}(object, ...)

}
\arguments{

  \item{...}{
    Additional arguments. Not used in any method.
  }

  \item{attrs}{
    Character vector of attribute names to be returned from the
    database. For a complete list of allowed attributes use the
    \code{\link{listColumns}} method.
  }

  \item{by}{
    For \code{mtisBy}: by which attribute should the MTIs be grouped?
    Allowed values are: \code{"gene"}, \code{"matmirna"},
    \code{"entrezid"}, \code{"pmid"}, \code{"support_type"},
    \code{"premirna"} and \code{"mirfam"}.
  }

  \item{filter}{
    A \code{list} of Filter classes (extending the \code{BasicFilter} class in
    the \code{ensembldb} package), like
    \code{\linkS4class{SpeciesFilter}} that should be used to filter the
    results returned from the database. See the documentation of the
    \code{\linkS4class{SpeciesFilter}} for a complete list of supported
    filters and the examples below for usage.
  }

  \item{force}{
    \code{force=TRUE} forces the function to return results, even if
    e.g. the pre-miRNAs specified in a \code{PremirnaFilter} could not
    be mapped to a mature miRNA name.
  }

  \item{object}{
    The \code{MirtarbaseDb} instance.
  }

  \item{of}{
    For \code{listSpecies}: either \code{"gene"} or
    \code{"mirna"} specifying whether species names for target genes or
    miRNAs should be returned.
  }

  \item{order.by}{
    On which attribute the results (if
    \code{return.type="data.frame"}) should be ordered. For
    \code{return.type="MTI"} the interactions will always be ordered
    by the ID of the MTI.
  }

  \item{order.type}{
    If the result (\code{data.frame}) should be ordered ascending
    (\code{order.type="asc"}) or descending (\code{order.type="desc"}).
  }

  \item{return.type}{
    Either \code{"data.frame"} or \code{"MTIList"} specifying whether the
    results should be returned as a \code{data.frame} or a \code{\linkS4class{MTIList}}
    (list of \code{\linkS4class{MTI}} objects).
  }

  \item{skip.keys}{
    For \code{listColumns}: whether primary and foreign keys should be
    returned.
  }

  \item{table}{
    For \code{listColumns}: specify the database table name for which
    the attributes should be listed.
  }

  \item{x}{
    The \code{MirtarbaseDb} instance.
  }

}
\section{Slots}{
  \describe{
    \item{\code{con}:}{
      Object of class \code{"DBIDbconn"}: the
      dbconn to the database.
    }

    \item{\code{mirtarbase_date}:}{
      The date when the miRTarbase database was created.
    }

    \item{\code{mirtarbase_version}:}{
      The version of the miRTarbase database.
    }

    \item{\code{species_mirna}:}{
      A character vector of all species names for miRNAs for which
      an MTI is defined in the database.
    }

    \item{\code{species_target_gene}:}{
      A character vector of all species names for target genes for which
      an MTI is defined in the database.
    }

    \item{\code{support_type}:}{
      A character vector of all support types of the MTIs.
    }

    \item{\code{tables}:}{
      named list of database table attributes with the names being the
      database table names. The order in which the tables are returned
      are the by their degree, i.e. the number of other tables they can
      be connected to (joined with).
    }

  }
}
\section{Methods and Functions}{
  \describe{

    \item{dbconn}{returns dbconn
      to the internal SQL database.
    }

    \item{listColumns}{
      Returns a character vector of all attributes stored in a database,
      or in a single database table if \code{table} was specified.
    }

    \item{listExperiments}{
      Returns a character vector of all experiments performed to
      validate the MTIs. Note that the naming of the experiments is
      quite inconsistent, e.g. real time qPCR is used as well as
      RT-qPCR.
    }

    \item{listPmids}{
      Returns a character vector of Pubmed identifiers for the
      publications in which the MTIs have been described.
    }

    \item{listSpecies}{
      Returns a character vector of the species names for either target
      genes or miRNAs for which an MTI is defined in the database.
    }

    \item{listSupportTypes}{
      Returns a character vector of all support types (weak and strong)
      of the MTIs.
    }

    \item{mtis}{
      Returns either a \code{\linkS4class{MTIList}} (list of \code{\linkS4class{MTI}}) objects or a
      \code{data.frame} listing all miRNA target gene interactions
      (MTIs) defined in the database passing the provided filters. See
      examples below for usage.
      Note that the returned \code{data.frame} will usually have more
      rows than actual MTIs, as each row in the \code{data.frame}
      represents one publication in which the MTI was validated. The
      returned \code{\linkS4class{MTIList}} on the other hand returns unique MTIs, as the
      supporting publications are stored in the \code{report} slot of
      the \code{\linkS4class{MTI}} objects.
    }

    \item{mtisBy}{
      Returns miRNA target gene interactions grouped by the attribute
      specified by \code{by}. The result is a \code{list} of
      \code{\linkS4class{MTIList}}s with the names of the list being the
      values of \code{by} and the elements the list of MTIs.
    }

    \item{show}{
      displays some informations from the database.
    }

    \item{listTables}{
      returns a named list of database table attributes (names of the
      list being the database table names).
    }

    \item{version}{
      returns the version of the miRTarbase resource
      that was used to build the database.
    }

  }
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{SpeciesFilter}}, \code{\linkS4class{MTI}}, \code{\linkS4class{MTIList}}
}
\examples{

## DB information
mirtarbase

## version
version(mirtarbase)

## what attributes do we have?
listColumns(mirtarbase)

## what tables
listTables(mirtarbase)

## list all species for target genes
listSpecies(mirtarbase, "gene")

## list all species for miRNAs
listSpecies(mirtarbase, "mirna")

## list all MTI support types
listSupportTypes(mirtarbase)

## list all experiments
listExperiments(mirtarbase)

## perform a SQL query on the database.
dbGetQuery(dbconn(mirtarbase), "select * from mirtarbase limit 2")


##******************************************
##
##   mtis:
##   retrieve miRNA target gene interactions
##   from the database.
##
##******************************************
## get all MTIs for gene BCL2 and restrict to human.
BCL2 <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2"), SpeciesFilter("Homo sapiens", feature="gene")))

## for how many miRNAs there is evidence that they target this gene?
length(BCL2)

BCL2
## just looking at the first
BCL2[[1]]

## what evidences are there for the interaction?
table(unlist(supportedBy(BCL2)))

## which pre-miRNAs are these?
pres <- premirna(BCL2)
BCL2[ which(is.na(pres)) ]
## well, some of the mature miRNAs can not be mapped to pre-miRNAs, because
## the mirbase version between the mirbase.db and the miRTarbase does not
## match.

## is there any "enrichment" for a specific miRNA family?
sort(table(mirfam(BCL2)), decreasing=TRUE)

## we could also use the lapply to apply the mirfam method to each
## individual MTI object.
sort(table(unlist(lapply(BCL2, mirfam))), decreasing=TRUE)

## now we perform the same query as above, but ask for a data.frame as
## return type.
BCL2.df <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2"), SpeciesFilter("Homo sapiens", feature="gene")), return.type="data.frame")

## the query is much faster, but:
nrow(BCL2.df)
length(BCL2)
## we get almost twice as many rows: one row in the data.frame corresponds to one
## publication in which a miRNA target gene interaction was detected, while each
## element in the list corresponds to one MTI class with the list of Report objects
## (i.e. publications) being stored in the report slot.

## get all MTIs for a specific miRNA
Test <- mtis(mirtarbase, filter=list(PremirnaFilter(c("hsa-mir-223")), SpeciesFilter("Homo sapiens")))
length(Test)
## get the gene names along with the number of supporting evidences
data.frame(gene=gene(Test), report_count=reportCount(Test))


## get all MTIs for a miRNA family.
MTIs <- mtis(mirtarbase, filter=list(MirfamFilter("mir-15"), SpeciesFilter("Homo sapiens")))
length(MTIs)
MTIs
## next we might want to ask if there are genes targeted by more than one miRNA
## of this family.
head(sort(table(gene(MTIs)), decreasing=TRUE))


}
\keyword{classes}

