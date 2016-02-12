## Package initialization:
.onLoad <- function(libname, pkgname){
    ns <- asNamespace(pkgname)
    path <- system.file("extdata", package=pkgname)
    objname <- "mirtarbase"
    db <- mirtarbaseDb(paste0(path, "/", objname, ".sqlite"))
    assign(objname, db, envir=ns)
    tmp <- getSpeciesDF()  ## call this once, so that we read that.
    assign(".mti.species", tmp, envir=ns)
    namespaceExport(ns, objname)
    ## setting useFancyQuotes to FALSE!
    .mti.fq <- options("useFancyQuotes")[[ 1 ]]
    options(useFancyQuotes=FALSE)
}

.onAttach <- function(libname, pkgname){
    packageStartupMessage(paste("Note: global option \"useFancyQuotes\" was set to FALSE.\n"))
}

.onUnload <- function(libpath){
    ## disconnecting from the database
    Vars <- ls(pos=1, all.names=TRUE)
    ## restoring the useFancyQuotes to the pre-load state.
    if(any(Vars == ".mti.fq")){
        options(useFancyQuotes=get(".mti.fq"))
    }
    if(any(Vars == ".mti.con")){
        con <- get(".mti.con")
        dbDisconnect(con)
    }
}

