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

## DB information
mirtarbase

## version
version( mirtarbase )

## what attributes do we have?
listAttributes( mirtarbase )

## what tables
tables( mirtarbase )

## list all species for target genes
listSpecies( mirtarbase, "target_gene" )

## list all species for miRNAs
listSpecies( mirtarbase, "mirna" )

## list all MTI support types
listSupportTypes( mirtarbase )

## list all experiments
listExperiments( mirtarbase )

## species data frame.
tmp <- mirtarbase:::getSpeciesDF()
head( tmp )

## perform a SQL query on the database.
dbGetQuery( connection( mirtarbase ), "select * from mirtarbase limit 2" )
