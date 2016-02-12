\name{premirna2matmirna}
\alias{premirna2matmirna}
\alias{premirna2matmirnaAcc}
\alias{premirna2premirnaAcc}
\alias{premirna2mirfam}
\alias{premirna2mirfamAcc}
\alias{premirnaAcc2matmirna}
\alias{premirnaAcc2matmirnaAcc}
\alias{premirnaAcc2premirna}
\alias{premirnaAcc2mirfam}
\alias{premirnaAcc2mirfamAcc}
\alias{matmirna2premirna}
\alias{matmirna2premirnaAcc}
\alias{matmirna2matmirnaAcc}
\alias{matmirna2mirfam}
\alias{matmirna2mirfamAcc}
\alias{matmirnaAcc2premirna}
\alias{matmirnaAcc2premirnaAcc}
\alias{matmirnaAcc2matmirna}
\alias{matmirnaAcc2mirfam}
\alias{matmirnaAcc2mirfamAcc}
\alias{mirfam2premirna}
\alias{mirfam2premirnaAcc}
\alias{mirfam2mirfamAcc}
\alias{mirfam2matmirna}
\alias{mirfam2matmirnaAcc}
\alias{mirfamAcc2premirna}
\alias{mirfamAcc2premirnaAcc}
\alias{mirfamAcc2mirfam}
\alias{mirfamAcc2matmirna}
\alias{mirfamAcc2matmirnaAcc}

\title{
  Conversion between various miRNA identifiers
}
\description{

  These functions use the \code{mirbase.db} package to perform
  conversions between miRNA identifiers (miRNA name/id to miRNA
  accession number) or to retrieve mature miRNAs for pre-miRNAs etc.

}
\usage{

matmirna2matmirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirna2mirfam( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirna2mirfamAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirna2premirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirna2premirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirnaAcc2matmirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirnaAcc2mirfam( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirnaAcc2mirfamAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirnaAcc2premirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
matmirnaAcc2premirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )

mirfam2matmirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfam2matmirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfam2mirfamAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfam2premirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfam2premirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfamAcc2matmirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfamAcc2matmirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfamAcc2mirfam( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfamAcc2premirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
mirfamAcc2premirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )

premirna2matmirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirna2matmirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirna2mirfam( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirna2mirfamAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirna2premirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirnaAcc2matmirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirnaAcc2matmirnaAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirnaAcc2mirfam( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirnaAcc2mirfamAcc( x, condition="=", ifnotfound=NA, return.type="data.frame" )
premirnaAcc2premirna( x, condition="=", ifnotfound=NA, return.type="data.frame" )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    a character vector or a single identifier.
  }

  \item{condition}{
    The condition to search for the values. Allowed are: \code{"="},
    \code{"!="}, \code{"like"} and \code{"not like"}, but in most cases
    the default (\code{"="}) is the most meaningful condition.
  }

  \item{ifnotfound}{
    Value to be returned if nothing was found for \code{x}.
  }

  \item{return.type}{
    Whether the returned object should be a \code{data.frame}
    (\code{return.type="data.frame"}) or a \code{list}
    (\code{return.type="list"}).
  }

}
\value{
  A \code{data.frame} with two columns, the first containing the query,
  the second the values, or a \code{list} of values with the names of
  the elements corresponding to the query.
}
\details{
  These function retrieve informations from the \code{mirbase.db}
  Bioconductor package.
}
\author{
Johannes Rainer
}
\examples{

## pre-miRNA to mature miRNA
Res <- premirna2matmirna( c( "bla", "hsa-mir-16-1", "hsa-mir-15b", "hsa-mir-16-2" ), return.type="list" )
Res
Res <- premirna2matmirnaAcc( c( "hsa-mir-15b" ) )
Res
Res <- premirnaAcc2matmirna( c( "MI0000070" ) )
Res
Res <- premirnaAcc2matmirnaAcc( c( "MI0000070" ) )
Res
Res <- premirna2premirnaAcc( "hsa-mir-15b" )
Res
Res <- premirnaAcc2premirna( "MI0000070" )
Res

## mature miRNA to pre-miRNA
Res <- matmirna2premirna( c( "hsa-miR-16-5p", "hsa-miR-16-3p" ) )
Res
Res <- matmirna2premirnaAcc( c( "hsa-miR-16-5p", "hsa-miR-16-3p" ) )
Res
Res <- matmirnaAcc2premirna( c( "MIMAT0000069", "MIMAT0004489" ) )
Res
Res <- matmirnaAcc2premirnaAcc( "MIMAT0004489" )
Res
Res <- matmirna2matmirnaAcc( "hsa-miR-16-5p" )
Res
Res <- matmirnaAcc2matmirna( "MIMAT0004489" )
Res

## pre-miRNA to mirfam and vice versa
Res <- premirna2mirfam( c( "hsa-mir-16-1", "hsa-mir-15b", "hsa-mir-16-2" ) )
Res
Res <- premirnaAcc2mirfam( "MI0000070" )
Res
Res <- mirfam2premirna( "mir-15" )
head( Res )
Res <- mirfam2premirnaAcc( "mir-15" )
head( Res )
Res <- mirfam2mirfamAcc( "mir-15" )
Res
Res <- mirfamAcc2mirfam( "MIPF0000006" )
Res
Res <- mirfamAcc2premirna( "MIPF0000006" )
head( Res )


## mature miRNAs to mirfam and vice versa
Res <- matmirna2mirfam( "hsa-miR-15b-5p" )
Res
Res <- matmirna2mirfamAcc( c( "bla", "hsa-miR-15b-5p" ) )
Res
Res <- mirfam2matmirna( "mir-15" )
head( Res )
Res <- mirfamAcc2matmirna( "MIPF0000006" )
head( Res )
Res <- mirfam2matmirnaAcc( "mir-15" )
head( Res )

}
\keyword{ data }
