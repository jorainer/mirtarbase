\name{ExperimentFilter}
\alias{ExperimentFilter}
\alias{PublicationFilter}
\alias{SpeciesFilter}
\alias{SupportTypeFilter}
\title{
  Constructor functions for filter objects
}
\description{
  These functions allow to create filter objects that can be used to
  retrieve specific elements from the annotation database.
}
\usage{
ExperimentFilter(value, condition = "=" )
PublicationFilter(value, condition = "=" )
SpeciesFilter(value, condition = "=", feature="gene" )
SupportTypeFilter(value, condition = "=" )
}
\arguments{
  \item{condition}{
    The condition to be used in the comparison. For character values
    \code{"="}, \code{"in"} and \code{"like"} are allowed, for numeric values
    (\code{SeqstartFilter} and \code{SeqendFilter}) \code{"="},
    \code{">"}, \code{">="}, \code{"<"} and \code{"<="}. Note that for
    \code{"like"} \code{value} should be a SQL pattern (e.g. \code{"ENS\%"}).
  }

  \item{feature}{
    For \code{SpeciesFilter}: allows to specify whether entries should
    be filtered by the genes of a certain species (\code{"gene"}) or
    miRNAs from a certain species (\code{"mirna"}).
  }

  \item{value}{
    The filter value, e.g., for \code{PublicationFilter} the Pubmed ID of the feature
    for which the data should be retrieved. Use functions
    \code{\link{listExperiments}}, \code{\link{listSpecies}} and
    \code{\link{listSupportTypes}} to retrieve all possible values for
    the respective filters.
  }
}
\value{
  Depending on the function called an instance of:
  \code{\linkS4class{ExperimentFilter}},
  \code{\linkS4class{PublicationFilter}},
  \code{\linkS4class{SpeciesFilter}},
  \code{\linkS4class{SupportTypeFilter}}
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{ExperimentFilter}},
  \code{\linkS4class{PublicationFilter}},
  \code{\linkS4class{SpeciesFilter}},
  \code{\linkS4class{SupportTypeFilter}}
  \code{\link[ensembldb]{GenenameFilter}}
}
\examples{

## for examples see examples for Filter classes

}
\keyword{data}
