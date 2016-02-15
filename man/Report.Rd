\name{MTIReport}
\alias{MTIReport}
\title{
  Create a Report object.
}
\description{
  This function creates a new \code{\linkS4class{MTIReport}} object.
}
\usage{
MTIReport(pmid=0, support_type="", experiments="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pmid}{the PubMed identifier of the report/paper.}
  \item{support_type}{the support type.}
  \item{experiments}{a character vector of experiments.}
}
\value{
  A MTIReport object.
}
\seealso{
  \code{\linkS4class{MTI}}, \code{\linkS4class{MTIReport}}
}
\author{
Johannes Rainer
}
\keyword{ data }

