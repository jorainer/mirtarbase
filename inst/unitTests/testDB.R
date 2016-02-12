## some basic tests for methods of the database.
## detachem <- function(x){
##     NS <- loadedNamespaces()
##     if(any(NS==x)){
##         pkgn <- paste0("package:", x)
##         detach(pkgn, unload=TRUE, character.only=TRUE)
##     }
## }
## Pkgs <- c("mirtarbase")
## tmp <- sapply(Pkgs, detachem)
## tmp <- sapply(Pkgs, library, character.only=TRUE)

test_loadDb <- function(){

    ## DB information
    mirtarbase

    ## version
    version(mirtarbase)

    ## what attributes do we have?
    listColumns(mirtarbase)

    ## what tables
    listTables(mirtarbase)

    ## list all species for target genes
    checkTrue(any(listSpecies(mirtarbase, "gene")=="Homo sapiens"))

    ## list all species for miRNAs
    listSpecies(mirtarbase, "mirna")

    ## list all MTI support types
    listSupportTypes(mirtarbase)

    ## list all experiments
    listExperiments(mirtarbase)

    ## species data frame.
    tmp <- mirtarbase:::getSpeciesDF()
    head(tmp)

    ## perform a SQL query on the database.
    dbGetQuery(dbconn(mirtarbase), "select * from mirtarbase limit 2")
}

