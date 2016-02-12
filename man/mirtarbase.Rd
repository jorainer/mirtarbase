\name{mirtarbase}
\alias{mirtarbase}
\title{
  miRTarBase (MTI) database package.
}
\description{
  This package contains the miRNA-target gene interactions (MTI)
  provided by the miRTarBase.
}
\usage{
  mirtarbase
}
\details{
  The \code{\linkS4class{MirtarbaseDb}} object representing the
  connection to the (internal) miRTarbase database.
}
\author{
Johannes Rainer
}
\seealso{
  \code{\link{mtis}}, \code{\link{listSupportTypes}}, \code{\link{listSpecies}}
}
\examples{

## print some summary information for the database:
mirtarbase

## what type of object did we get?
class( mirtarbase )

## list all available evidence support types from the database
listSupportTypes( mirtarbase )

## list all species (for target genes):
listSpecies( mirtarbase )

}
\keyword{ data }

