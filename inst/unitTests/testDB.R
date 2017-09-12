## some basic tests for methods of the database.

test_loadDb <- function(){
    ## DB information
    mirtarbase
    ## version
    checkTrue(length(version(mirtarbase)) > 0)
    ## what attributes do we have?
    checkTrue(length(listColumns(mirtarbase)) > 0)
    ## what tables
    checkTrue(length(listTables(mirtarbase)) > 0)
    ## list all species for target genes
    checkTrue(any(listSpecies(mirtarbase, "gene")=="Homo sapiens"))
    ## list all species for miRNAs
    checkTrue(length(listSpecies(mirtarbase, "mirna")) > 0)
    ## list all MTI support types
    checkTrue(length(listSupportTypes(mirtarbase)) > 0)
    ## list all experiments
    checkTrue(length(listExperiments(mirtarbase)) > 0)
    ## species data frame.
    tmp <- mirtarbase:::getSpeciesDF()
    head(tmp)
    checkTrue(nrow(tmp) > 0)    
}

