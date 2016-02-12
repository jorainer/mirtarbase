\name{PositionFilter-class}
\Rdversion{1.1}
\docType{class}
\alias{ExperimentFilter-class}
\alias{PublicationFilter-class}
\alias{SpeciesFilter-class}
\alias{SupportTypeFilter-class}
% defined in ensembldb
\alias{attribute,EntrezidFilter,MirtarbaseDb-method}
\alias{where,EntrezidFilter,MirtarbaseDb-method}
\alias{attribute,GenenameFilter,MirtarbaseDb-method}
\alias{where,GenenameFilter,MirtarbaseDb-method}
% defined in mirnahostgenes
\alias{attribute,MatmirnaFilter,MirtarbaseDb-method}
\alias{where,MatmirnaFilter,MirtarbaseDb-method}
\alias{attribute,MatmirnaidFilter,MirtarbaseDb-method}
\alias{where,MatmirnaidFilter,MirtarbaseDb-method}
\alias{attribute,MirfamFilter,MirtarbaseDb-method}
\alias{where,MirfamFilter,MirtarbaseDb-method}
\alias{attribute,MirfamidFilter,MirtarbaseDb-method}
\alias{where,MirfamidFilter,MirtarbaseDb-method}
\alias{attribute,PremirnaFilter,MirtarbaseDb-method}
\alias{where,PremirnaFilter,MirtarbaseDb-method}
\alias{attribute,PremirnaidFilter,MirtarbaseDb-method}
\alias{where,PremirnaidFilter,MirtarbaseDb-method}
% defined in mirtarbase
\alias{attribute,ExperimentFilter,MirtarbaseDb-method}
\alias{attribute,PublicationFilter,MirtarbaseDb-method}
\alias{attribute,SpeciesFilter,MirtarbaseDb-method}
\alias{attribute,SupportTypeFilter,MirtarbaseDb-method}
\alias{where,ExperimentFilter,MirtarbaseDb-method}
\alias{where,PublicationFilter,MirtarbaseDb-method}
\alias{show,SpeciesFilter,MirtarbaseDb-method}
\alias{show,SpeciesFilter-method}
\alias{where,SpeciesFilter,MirtarbaseDb-method}
\alias{where,SupportTypeFilter,MirtarbaseDb-method}
\alias{attribute}
\alias{where}

\title{Filter results fetched from the miRTarbase database}
\description{
  These classes allow to specify which entries should be retrieved from the database.
}
\section{Objects from the Class}{
  While objects can be created by calls e.g. of the form
  \code{new("ExperimentFilter", ...)} users are strongly encouraged to use
  the specific functions: \code{\link{ExperimentFilter}},
  \code{\link{PublicationFilter}}, \code{\link{SpeciesFilter}} and
  \code{\link{SupportTypeFilter}}.
  See examples below for usage.
}
\section{Slots}{
  \describe{
    \item{\code{condition}:}{
      Object of class \code{"character"}: can be
      either \code{"="}, \code{"in"} or \code{"like"} to filter on character values
      (e.g. gene id, gene biotype, seqname etc), or \code{"="}, \code{">"}
      or \code{"<"} for numerical values (chromosome/seq
      coordinates). Note that for \code{"like"} \code{value} should be a
      SQL pattern (e.g. \code{"ENS\%"}).
    }

    \item{feature}{
      For \code{SpeciesFilter}: allows to specify whether entries should
      be filtered by the genes of a certain species (\code{"gene"}) or
      miRNAs from a certain species (\code{"mirna"}).
    }

    \item{\code{value}:}{
      Object of class \code{"character"}: the value
      to be used for filtering. Possible values for
      \code{ExperimentFilter}, \code{SpeciesFilter} or
      \code{SupportTypeFilter} can be retrieved using the
      \code{\link{listExperiments}}, \code{\link{listSpecies}} and
      \code{\link{listSupportTypes}}, respectively.
    }

    \item{\code{valueIsCharacter}:}{
      Object of class \code{"logical"}:
      that's an internal slot which value should not be changed.
    }
  }
}
\section{Extends}{
Class \code{"BasicFilter"} defined in the \code{ensembldb} package, directly.
}
\section{Methods}{
  \describe{

    \item{attribute}{
      \code{signature(object = "ExperimentFilter", db = "MirtarbaseDb")}: returns
      the attribute (table column name in the miRTarbase database)
      to be used for the filtering. The method can be applied to any
      filter class defined in the \code{mirtarbase}, \code{mirnahostgenes} or
      \code{ensembldb} packages.
    }

    \code{show}{
      Print some information about the object.
    }

    \item{where}{
      \code{signature(object = "ExperimentFilter", db = "MirtarbaseDb")}: returns the
      where clause for the SQL call to the miRTarbase database. The method can be applied to any
      filter class defined in the \code{mirtarbase}, \code{mirnahostgenes} or
      \code{ensembldb} packages.
    }

  }
}
\details{
  In addition to the above classes all Filter classes defined in the
  \code{ensembldb} or \code{mirnahostgenes} packages.
  The complete list of allowed filters is (in alphabetic order):
  \describe{

    \item{\code{EntrezidFilter}}{
      Filter the results for specific NCBI Entrezgene IDs.
    }

    \item{\code{ExperimentFilter}}{
      Restrict the results to specific experiments. See
      \code{\link{listExperiments}} for a full list of supported
      experiments.
    }

    \item{\code{GenenameFilter}}{
      Restrict the results to specific genes based on their name/symbol.
    }

    \item{\code{MatmirnaFilter}}{
      Restrict the results to specified mature miRNA names.
    }

    \item{\code{MatmirnaidFilter}}{
      Restrict the results to specified mature miRNA accession IDs.
    }

    \item{\code{MirfamFilter}}{
      Restrict the results to interactions of miRNAs from the specified
      miRNA family.
    }

    \item{\code{MirfamidFilter}}{
      Restrict the results to interactions of miRNAs from the specified
      miRNA family (ID/accession).
    }

    \item{\code{PremirnaFilter}}{
      Restrict the results to interactions of mature miRNAs derived from
      the specified pre-miRNAs (identified by their name).
    }

    \item{\code{PremirnaidFilter}}{
      Restrict the results to interactions of mature miRNAs derived from
      the specified pre-miRNAs (identified by their accession).
    }

    \item{\code{PublicationFilter}}{
      Restrict the results to MTIs validated in the specified
      publications (identified by their PubMed IDs).
    }

    \item{\code{SpeciesFilter}}{
      Restrict the results to miRNAs or target genes from the specified species.
    }

    \item{\code{SupportTypeFilter}}{
      Restrict the results to MTIs with specified support type/evidence
      level. Use the method \code{\link{listSupportTypes}} to get a list
      of available support types.
    }

  }
}
\note{
  The filter for experiments should ideally use as condition a
  \code{"like"} and as value a pattern such as \code{"\%PCR\%"} to retrieve
  all entries for MTIs validated using PCR based methods as the
  values in the experiments column of the miRTarbase database are not
  standardized.

  The miRTarbase database does only provide mapping between mature miRNA
  names and target genes, thus, \code{PremirnaFilter},
  \code{PremirnaidFilter}, \code{MatmirnaidFilter}, \code{MirfamFilter}
  and \code{MirfamidFilter} instances internally map the provided
  identifier to the corresponding mature miRNA names using the
  \code{mirbase.db} package (e.g. using the
  \code{\link{premirna2matmirna}} function).
}
\author{
  Johannes Rainer
}
\seealso{
  \code{\link[ensembldb]{GenenameFilter}}, \code{\link{premirna2matmirna}}
}
\examples{

##***********************
##
## EntrezidFilter
##***********************
##
EFilt <- EntrezidFilter( 123, condition="like" )
EFilt
## default attribute of the filter
attribute( EFilt )
## attribute of the filter for the mirtarbase package:
attribute( EFilt, mirtarbase )
## and the where query:
where( EFilt, mirtarbase )

## EntrezidFilter with more than one value:
EFilt <- EntrezidFilter( c( 123, 3435 ) )
EFilt
where( EFilt, mirtarbase )


##***********************
##
## ExperimentFilter
##***********************
## create a simple experiment filter; this would return all entries that have "qpcr"
## in the experiments attribute; like is case insensitive.
EFilt <- ExperimentFilter( "PCR" )
attribute( EFilt, mirtarbase )
where( EFilt, mirtarbase )

EFilt <- ExperimentFilter( c( "PCR" ), condition="like" )
where( EFilt, mirtarbase )

## a good starting point is to list all experiments from the database
## and select the pattern from there
listExperiments( mirtarbase )

##***********************
##
## GenenameFilter
##***********************
## create a gene name filter
GFilt <- GenenameFilter( "Bcl2l11", condition="!=" )
GFilt
## the attribute name in the ensembldb package
attribute( GFilt )
## and in the mirtarbase package
attribute( GFilt, mirtarbase )
where( GFilt, mirtarbase )


##***********************
##
## MatmirnaFilter
##***********************
## mature miRNA filter
MFilt <- MatmirnaFilter( c( "hsa-miR-16-5p", "hsa-miR-16-3p" ) )
MFilt
## and in the mirtarbase package
attribute( MFilt, mirtarbase )
where( MFilt, mirtarbase )


##***********************
##
## PublicationFilter
##***********************
## create a simple publication filter (on Pubmed IDs)
PFilt <- PublicationFilter( 123 )
PFilt
## the attribute in the mirtarbase package
attribute( PFilt, mirtarbase )
where( PFilt, mirtarbase )


##***********************
##
## SpeciesFilter
##***********************
## create a species filter on genes
SFilt <- SpeciesFilter( "Homo sapiens", feature="gene" )
SFilt
## the attribute for the gene
attribute( SFilt, mirtarbase )
where( SFilt, mirtarbase )

## the same for miRNA
SFilt <- SpeciesFilter( "Homo sapiens", feature="mirna" )
## the attribute for miRNA
attribute( SFilt, mirtarbase )
where( SFilt, mirtarbase )

## what if we used a non-standard species?
SFilt <- SpeciesFilter( "Homo hobbeniensis", feature="gene" )
## this throws a warning
suppressWarnings(
    where( SFilt, mirtarbase )
)
## obviously only values present in the database make sense:
listSpecies( mirtarbase, "gene" )
listSpecies( mirtarbase, "mirna" )


##***********************
##
## SupportTypeFilter
##***********************
## create a filter on the mirtarbase support types
SFilt <- SupportTypeFilter( "unknown" )
SFilt
## the attribute that will be queried
attribute( SFilt, mirtarbase )
## throws a warning since the support type unknown is not present in the database
suppressWarnings(
    where( SFilt, mirtarbase )
    )

## better to select one of the correct support types:
listSupportTypes( mirtarbase )


##***********************
##
## PremirnaFilter
##***********************
## pre miRNA filter.
## Note that the PremirnaFilter methods for MirtarbaseDb internally map the pre-miRNA name
## to mature miRNA names based on the mirbase.db package, as the
## miRTarbase database does not provide pre-miRNA ids.
PF <- PremirnaFilter( c( "bla", "hsa-mir-15b" ) )
## the where internally maps the pre-miRNA ID to the corresponding mature miRNA ID(s)
where( PF, mirtarbase )

## what if nothing was found...
PF <- PremirnaFilter( "baf" )
where( PF, mirtarbase )
## returns nothing

PF <- PremirnaFilter( c( "bla", "hsa-mir-15b" ), condition="!=" )
where( PF, mirtarbase )

##
PF <- PremirnaFilter( "hsa-mir-15b" )
where(PF , mirtarbase)


##***********************
##
## PremirnaidFilter
##***********************
PF <- PremirnaidFilter( "MI0000070" )
where( PF, mirtarbase )


##***********************
##
## MatmirnaidFilter
##***********************
MF <- MatmirnaidFilter( "MIMAT0004489", "!=" )
where( MF, mirtarbase )


##***********************
##
## MirfamFilter
##***********************
MF <- MirfamFilter( "mir-15" )
where( MF, mirtarbase )


##***********************
##
## MirfamidFilter
##***********************
MF <- MirfamidFilter( "MIPF0000006" )
where( MF, mirtarbase )


}
\keyword{classes}

