## Package initialization:
.onLoad <- function(libname, pkgname){
    ns <- asNamespace(pkgname)
    path <- system.file("extdata", package=pkgname)
    files <- dir(path, pattern="sqlite")
    for(i in seq_len(length(files))){
        db <- MirtarbaseDb(system.file("extdata", files[[i]], package=pkgname,
                                       lib.loc=libname))
        objname <- sub(".sqlite$","",files[[i]])
        assign(objname, db, envir=ns)
        namespaceExport(ns, objname)
    }
    versions <- gsub(files, pattern="MirtarbaseDb.v", replacement="", fixed=TRUE)
    versions <- gsub(versions, pattern=".sqlite", replacement="", fixed=TRUE)
    versions <- sort(versions, decreasing=TRUE)
    objname <- "mirtarbase"
    assign(objname, get(paste0("MirtarbaseDb.v", versions[1]), envir=ns), envir=ns)
    tmp <- getSpeciesDF()  ## call this once, so that we read that.
    assign(".mti.species", tmp, envir=ns)
    namespaceExport(ns, objname)
    ## setting useFancyQuotes to FALSE!
    .mti.fq <- options("useFancyQuotes")[[ 1 ]]
    options(useFancyQuotes=FALSE)
}

.onAttach <- function(libname, pkgname){
    ## Loaded xx mirtarbase versions. The shortcut mirtarbase links to the most recent one (6.1).
    path <- system.file("extdata", package=pkgname)
    files <- dir(path, pattern="sqlite")
    versions <- gsub(files, pattern="MirtarbaseDb.v", replacement="")
    versions <- gsub(versions, pattern=".sqlite", replacement="", fixed=TRUE)
    versions <- sort(versions, decreasing=TRUE)
    packageStartupMessage(paste0("mirtarbase:\nLoaded ", length(versions), " miRTarBase release(s).\n",
                                 "The shortcut 'mirtarbase' links to the most recent one (", versions[1],")\n",
                                 "Note: global option \"useFancyQuotes\" was set to FALSE.\n"))
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

