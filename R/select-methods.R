####============================================================
##  Methods and functions for/from AnnotationDbi
##
##
####------------------------------------------------------------

####============================================================
##  columns method.
##
####------------------------------------------------------------
setMethod("columns", "MirtarbaseDb",
          function(x).getColumns(x))
.getColumns <- function(x){
    return(names(.MIRTARBASE_COLMAP))
}
## Mapping vector: names are "AnnotationDbi-conform" uppercase names, elements
## are the columns from the mirtarbase database.
.MIRTARBASE_COLMAP <- c("mirtarbase_id", "mirna", "species_mirna", "target_gene",
                        "target_gene_entrez_gene_id", "species_target_gene",
                        "experiments", "support_type", "references_pmid")
names(.MIRTARBASE_COLMAP) <- c("MIRTARBASEID", "MATMIRNA", "MIRNASPECIES",
                               "SYMBOL", "ENTREZID", "GENESPECIES", "EXPERIMENT",
                               "SUPPORTTYPE", "PMID")
## Map AnnotationDbi formatted columns to MirtarbaseDb column names.
.mapCols2Mirtarbase <- function(x){
    if(!all(x %in% names(.MIRTARBASE_COLMAP))){
        whichNotAllowed <- x[!(x %in% names(.MIRTARBASE_COLMAP))]
        stop("Columns ", paste0(whichNotAllowed, collapse=", "), " are not known!",
             " Use the 'columns' method to list all supported column names.")
    }
    return(unname(.MIRTARBASE_COLMAP[x]))
}
## Other way round
.mapMirtarbase2Cols <- function(x){
    if(!all(x %in% .MIRTARBASE_COLMAP)){
        whichNotAllowed <- x[!(x %in% .MIRTARBASE_COLMAP)]
        stop("Don't know how to map ", paste0(whichNotAllowed, collapse=", "), "!")
    }
    tmp <- names(.MIRTARBASE_COLMAP)
    names(tmp) <- .MIRTARBASE_COLMAP
    return(unname(tmp[x]))
}

####============================================================
##  keytypes.
##
##  I'll use all of the supported Filters here.
##  EntrezidFilter
##  GenenameFilter
##  MatmirnaFilter
##  SpeciesFilter
##  ExperimentFilter
##  PublicationFilter
##  SupportTypeFilter
##  I'll skip the MatmirnaidFilter, PremirnaFilter, PremirnaidFilter,
##  MirfamFilter and MirfamidFilter for now.
####------------------------------------------------------------
setMethod("keytypes", "MirtarbaseDb",
          function(x){
              return(names(.keytype2FilterMapping()))
          }
)
## returns a vector mapping keytypes (names of vector) to filter names (elements).
.keytype2FilterMapping <- function(){
    filters <- c("EntrezidFilter", "GenenameFilter", "MatmirnaFilter", "SpeciesFilter",
                 "SpeciesFilter", "ExperimentFilter", "PublicationFilter", "SupportTypeFilter",
                 "MirtarbaseidFilter")
    names(filters) <- c("ENTREZID", "SYMBOL", "MATMIRNA", "GENESPECIES", "MIRNASPECIES",
                        "EXPERIMENT", "PMID", "SUPPORTTYPE", "MIRTARBASEID")
    return(filters)
}
filterForKeytype <- function(keytype){
    if(length(keytype) > 1){
        keytype <- keytype[1]
        warning("Multiple keytype values are not supported! Using only the first keytype.")
    }
    filters <- .keytype2FilterMapping()
    if(any(names(filters) == keytype)){
        ## Handle "special" cases.
        if(keytype == "MIRNASPECIES")
            return(SpeciesFilter(value="tobeset", feature="mirna"))
        if(keytype == "GENESPECIES")
            return(SpeciesFilter(value="tobeset", feature="gene"))
        filt <- new(filters[keytype])
        return(filt)
    }else{
        stop("No filter for that keytype!")
    }
}

####============================================================
##  keys method
##
##  The keys method returns all of the keys for a specified keytype.
##  If keytype is not specified we're returning the mirtarbase ids.
####------------------------------------------------------------
setMethod("keys", "MirtarbaseDb",
          function(x, keytype, filter, ...){
              if(missing(keytype))
                  keytype <- "MIRTARBASEID"
              filter <- .checkFilter(filter)
              keyt <- keytypes(x)
              keytype <- match.arg(keytype, keyt)
              ## Map the keytype to the appropriate column name.
              dbCol <- .mapCols2Mirtarbase(keytype)
              ## Perform the query.
              res <- .getWhat(x, columns=dbCol, filter=filter)[, dbCol]
              return(unname(res))
          })

####============================================================
##  select method
##
##
####------------------------------------------------------------
setMethod("select", "MirtarbaseDb",
          function(x, keys, columns, keytype, ...){
              if(missing(keys))
                  keys <- NULL
              if(missing(columns))
                  columns <- NULL
              if(missing(keytype))
                  keytype <- NULL
              return(.select(x=x, keys=keys, columns=columns, keytype=keytype, ...))
          })
.select <- function(x, keys=NULL, columns=NULL, keytype=NULL, ...){
    extraArgs <- list(...)
    ## Perform argument checking:
    ## columns:
    if(missing(columns) | is.null(columns))
        columns <- columns(x)
    notAvailable <- !(columns %in% columns(x))
    if(all(notAvailable))
        stop("None of the specified columns are avaliable in the database!")
    if(any(notAvailable)){
        warning("The following columns are not available in the database and have",
                " thus been removed: ", paste(columns[notAvailable], collapse=", "))
        columns <- columns[!notAvailable]
    }
    ## keys:
    if(is.null(keys) | missing(keys)){
        ## Get everything from the database...
        keys <- list()
    }else{
        if(!(is(keys, "character") | is(keys, "list") | is(keys, "BasicFilter")))
            stop("Argument keys should be a character vector, an object extending BasicFilter ",
                 "or a list of objects extending BasicFilter.")
        if(is(keys, "list")){
            if(!all(vapply(keys, is, logical(1L), "BasicFilter")))
                stop("If keys is a list it should be a list of objects extending BasicFilter!")
        }
        if(is(keys, "BasicFilter")){
            keys <- list(keys)
        }
        if(is(keys, "character")){
            if(is.null(keytype)){
                stop("Argument keytype is mandatory if keys is a character vector!")
            }
            ## Check also keytype:
            if(!(keytype %in% keytypes(x)))
                stop("keytype ", keytype, " not available in the database.",
                     " Use keytypes method to list all available keytypes.")
            ## Generate a filter object for the filters.
            keyFilter <- filterForKeytype(keytype)
            value(keyFilter) <- keys
            keys <- list(keyFilter)
            ## Add also the keytype itself to the columns.
            if(!any(columns == keytype))
                columns <- c(keytype, columns)
        }
    }
    ## Map the columns to column names we have in the database.
    mtbCols <- .mapCols2Mirtarbase(columns)
    ## OK, now perform the query given the filters we've got.
    res <- .getWhat(x, columns=mtbCols, filter=keys)
    colnames(res) <- .mapMirtarbase2Cols(colnames(res))
    return(res[, columns])
}

