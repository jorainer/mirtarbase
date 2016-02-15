\name{mirtarbase}
\alias{mirtarbase}
\alias{MirtarbaseDb.v4.5}
\alias{listMirtarbaseReleases}
\title{
  miRTarBase (MTI) database package.
}
\description{
  This package contains the miRNA-target gene interactions (MTI)
  provided by the miRTarBase.
}
\usage{
  mirtarbase

  listMirtarbaseReleases()
}
\details{
  \describe{
    \item{mirtarbase}{
      The \code{\linkS4class{MirtarbaseDb}} object representing the
      connection to the (internal) miRTarbase database. The
      \code{mirtarbase} always links to the most recent miRTarbase release
      provided by the package (use \code{\link{listMirtarbaseReleases}}
      function to list all provided releases).
    }

    \item{listMirtarbaseReleases}{
      Function to list all miRTarbase releases provided by the
      package. Returns a \code{data.frame} listing the variable name and
      the corresponding version. The \code{mirtarbase} variable always
      links to the most recent release.
    }
  }
}
\author{
Johannes Rainer
}
\seealso{
  \code{\link{mtis}}, \code{\link{listSupportTypes}},
  \code{\link{listSpecies}}
}
\examples{

## Print some summary information for the database:
mirtarbase

## List all releases that are provided.
listMirtarbaseReleases()

## what type of object did we get?
class(mirtarbase)

## list all available evidence support types from the database
listSupportTypes(mirtarbase)

## list all species (for target genes):
listSpecies(mirtarbase)

}
\keyword{ data }

