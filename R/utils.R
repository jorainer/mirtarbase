
## this function checks if we've got the .MTB.SPECIES available, and if not it reads it from
## extdata.
getSpeciesDF <- function(){
    ns <- asNamespace("mirtarbase")
    if(any(ls(ns, all.names=TRUE) ==".mti.species"))
        return(get(".mti.species", envir=ns))
    DF <- read.table(system.file("extdata/organisms.txt.gz", package="mirtarbase"),
                     sep="\t", as.is=TRUE, header=TRUE, comment.char="")
    assign(".mti.species", DF, envir=ns)
    return(DF)
}

