## some basic tests for methods of the database.
detachem <- function( x ){
    NS <- loadedNamespaces()
    if( any( NS==x ) ){
        pkgn <- paste0( "package:", x )
        detach( pkgn, unload=TRUE, character.only=TRUE )
    }
}
Pkgs <- c( "mirtarbase" )
tmp <- sapply( Pkgs, detachem )
tmp <- sapply( Pkgs, library, character.only=TRUE )

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
EFilt <- ExperimentFilter( "%qPCR%", condition="like" )
attribute( EFilt, mirtarbase )
where( EFilt, mirtarbase )
## like with more than two values; not that like does not support multiple values,
## thus it is reverted to "in"
EFilt <- ExperimentFilter( c( "qPCR", "microarray" ), condition="like" )
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
## note that the PremirnaFilter methods for MirtarbaseDb internally map the pre-miRNA name
## to mature miRNA names.
PF <- PremirnaFilter( c( "bla", "hsa-mir-15b" ) )
where( PF )
## doesn't make much sense...
attribute( PF, mirtarbase )  ## that's the database column in which we're going to search
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


MF <- MatmirnaidFilter( c( "MIMAT0004489", "dfd" ), "in" )
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





