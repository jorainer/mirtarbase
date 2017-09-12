\name{select}
\Rdversion{1.1}
\alias{select}
\alias{select,MirtarbaseDb-method}
\alias{columns,MirtarbaseDb-method}
\alias{keys,MirtarbaseDb-method}
\alias{keytypes,MirtarbaseDb-method}
%\alias{mapIds,MirtarbaseDb-method}

\title{Integration into the AnnotationDbi framework}

\description{
  Several of the methods available for \code{AnnotationDbi} objects are
  also implemented for \code{MirtarbaseDb} objects. This enables to extract
  data from \code{MirtarbaseDb} objects in a similar fashion than from objects
  inheriting from the base annotation package class
  \code{AnnotationDbi}.
  In addition to the \emph{standard} usage, the \code{select} for
  \code{MirtarbaseDb} objects support also the filter framework of the
  mirtarbasedb (and ensembdb) package and thus allow to perform more
  fine-grained queries to retrieve data.
}
\usage{

\S4method{columns}{MirtarbaseDb}(x)
\S4method{keys}{MirtarbaseDb}(x, keytype, filter,...)
\S4method{keytypes}{MirtarbaseDb}(x)
%\S4method{mapIds}{MirtarbaseDb}(x, keys, column, keytype, ..., multiVals)
\S4method{select}{MirtarbaseDb}(x, keys, columns, keytype, ...)

}
\arguments{

  (In alphabetic order)

  %% \item{column}{
  %%   For \code{mapIds}: the column to search on, i.e. from which values
  %%   should be retrieved.
  %% }

  \item{columns}{
    For \code{select}: the columns from which values should be
    retrieved. Use the \code{columns} method to list all possible
    columns.
  }

  \item{keys}{
    The keys/ids for which data should be retrieved from the
    database. This can be either a character vector of keys/IDs, a
    single filter object extending \code{\linkS4class{AnnotationFilter}} or a
    list of such objects. For the latter case it is possible to combine
    filters and thus perform more specific queries.
  }

  \item{keytype}{
    For \code{select}: the type (column) that matches
    the provided keys. This argument does not have to be specified if
    argument \code{keys} is a filter object extending
    \code{\linkS4class{AnnotationFilter}} or a \code{list} of such objects.

    For \code{keys}: which keys should be returned from the database.
  }

  \item{filter}{
    For \code{keys}: either a single object extending
    \code{\linkS4class{AnnotationFilter}} or a list of such object to
    retrieve only specific keys from the database.
  }

  %% \item{multiVals}{
  %%   What should \code{mapIds} do when there are multiple values that
  %%   could be returned? Options are: \code{"first"}, \code{"list"},
  %%   \code{"filter"}, \code{"asNA"}. See
  %%   \code{\link[AnnotationDbi]{mapIds}} for a detailed description.
  %% }

  \item{x}{
    The \code{MirtarbaseDb} object.
  }

  \item{...}{
    Not used.
  }

}
\section{Methods and Functions}{
  \describe{

    \item{columns}{
      List all the columns that can be retrieved by the \code{select}
      methods. Note that these column names are different from the ones
      supported by the \code{\link{mtis}} method. These can be listed by the
      \code{\link{listColumns}} method.

      Returns a character vector of supported column names.
    }

    \item{keys}{
      Retrieves all keys from the column name specified with
      \code{keytype}. By default (if \code{keytype} is not provided) it
      returns all Mirtarbase IDs.

      Returns a character vector of IDs.
    }

    \item{keytypes}{
      List all supported key types (column names).

      Returns a character vector of key types.
    }

    %% \item{mapIds}{
    %%   Retrieve the mapped ids for a set of keys that are of a particular
    %%   keytype. Argument \code{keys} can be either a character vector of
    %%   keys/IDs, a single filter object extending
    %%   \code{\linkS4class{BasicFilter}} or a list of such objects. For
    %%   the latter, the argument \code{keytype} does not have to be
    %%   specified. Importantly however, if the filtering system is used,
    %%   the ordering of the results might not represent the ordering of
    %%   the keys.

    %%   The method usually returns a named character vector or, depending
    %%   on the argument \code{multiVals} a named list, with names
    %%   corresponding to the keys (same ordering is only guaranteed if
    %%   \code{keys} is a character vector).
    %% }

    \item{select}{
      Retrieve the data as a \code{data.frame} based on parameters for
      selected \code{keys}, \code{columns} and \code{keytype}
      arguments. Multiple matches of the keys are returned in one row
      for each possible match. Argument \code{keys} can be either a
      character vector of keys/IDs, a single filter object extending
      \code{\linkS4class{AnnotationFilter}} or a list of such objects. For
      the latter, the argument \code{keytype} does not have to be
      specified.

      Returns a \code{data.frame} with the column names corresponding to
      the argument \code{columns} and rows with all data matching the
      criteria specified with \code{keys}.
    }

  }
}

\value{
  See method description above.
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\linkS4class{MirtarbaseIdFilter}}
  \code{\link{listColumns}}
  \code{\link{mtis}}
}
\examples{

library(mirtarbase)

## List all supported keytypes.
keytypes(mirtarbase)

## List all supported columns for the select method.
columns(mirtarbase)

## List /real/ database column names.
listColumns(mirtarbase)

## Retrieve all keys corresponding to Mirtarbase IDs.
mtiids <- keys(mirtarbase, keytype="MIRTARBASEID")
length(mtiids)
head(mtiids)

## Retrieve all keys corresponding to gene names/symbol for which we do have
## a MTI of support type "Functional MTI"
syms <- keys(mirtarbase, keytype="SYMBOL", filter=SupportTypeFilter("Functional MTI"))
length(syms)
head(syms)

## select:
## Retrieve all MTIs for the selected genes and support type.
select(mirtarbase, keys=list(GenenameFilter(c("BCL2", "BCL2L11")),
                             SupportTypeFilter("Functional MTI")),
       columns=c("SYMBOL", "MATMIRNA", "PMID"))


}
\keyword{classes}





