##*************************************************************************
##
##     Public functions.
##
##*************************************************************************


##*************************************************************************
##
##     Internal functions.
##
##*************************************************************************
## create a connection to the database and return the MirtarbaseDb object.
mirtarbaseDb <- function(x){
    options(useFancyQuotes=FALSE)
    lite <- dbDriver("SQLite")
    con <- dbConnect(lite, dbname=x, flags=SQLITE_RO)
    tables <- dbListTables(con)
    ## read the attributes for these tables.
    Tables <- vector(length=length(tables), "list")
    for(i in 1:length(Tables)){
        Tables[[ i ]] <- colnames(dbGetQuery(con, paste0("select * from ",
                                                         tables[ i ], " limit 1")))
    }
    names(Tables) <- tables
    ## read the info file.
    info <- read.table(system.file("extdata/txt/INFO", package="mirtarbase"),
                       header=TRUE, as.is=TRUE, sep="\t")
    ## get also some additional info...
    species_tg <- dbGetQuery(con, "select distinct species_target_gene from mirtarbase;")[ , 1 ]
    species_mirna <- dbGetQuery(con, "select distinct species_mirna from mirtarbase;")[ , 1 ]
    st <- dbGetQuery(con, "select distinct support_type from mirtarbase")[ , 1 ]
    MDB <- new("MirtarbaseDb",
               con=con,
               tables=Tables,
               mirtarbase_version=info[ info$key=="release_version", "value" ],
               mirtarbase_date=info[ info$key=="release_date", "value" ],
               species_target_gene=species_tg,
               species_mirna=species_mirna,
               support_type=st
               )
    return(MDB)
}


## builds the query that we use to retrieve the data.
## x is the MirtarbaseDb object.
## filter is a list of filters.
## attrs: attributes to retrieve
## order.by
.buildQuery <- function(x, attrs=listAttributes(x), filter, order.by="",
                        order.type="asc", match.case=FALSE, force=FALSE){
    resultattrs <- attrs
    collatequery <- ""
    if(!missing(filter)){
        ## check filter!
        if(class(filter)!="list")
            stop("parameter filter has to be a list of BasicFilter classes!")
        ## add the attributes needed for the filter
        filterattrs <- unlist(lapply(filter, attribute, x))
        ##filterattrs <- sapply(filterattrs, removePrefix, USE.NAMES=FALSE)
        attrs <- unique(c(attrs, filterattrs))
        ## next we're building the where query.
        wheres <- lapply(filter, where, x)
        ## what if we've got any filter returning NULL in where?
        if(any(unlist(lapply(wheres, is.null))) & !force){
            ## this means we've got one of the mappings from e.g. pre-miRNA to mature
            ## miRNA completely failing. In such a case we would eventually fetch the
            ## complete mirtarbase.
            stop("One of the filters returned NULL. This most likely means that one of the submitted miRNA identifiers can not be mapped to a mature miRNA name! To still perform the query set parameter \"force=TRUE\".", call.=FALSE)
        }
        wherequery <- paste(" where", paste(unlist(wheres), collapse=" and "))
        ## should the query be performed case insensitive?
        if(!match.case){
            wherequery <-  paste0(wherequery, " collate nocase")
        }
    }else{
        wherequery <- ""
    }
    ## should we do an order.by?
    if(!missing(order.by) & order.by!=""){
        order.by <- unlist(strsplit(order.by, split=",", fixed=TRUE))
        order.by <- gsub(order.by, pattern=" ", replacement="", fixed=TRUE)
        ## allow only order.by that are also in the attributes.
        order.by.noattrs <- order.by[ !(order.by %in% attrs) ]
        order.by <- order.by[ order.by %in% attrs ]
        if(length(order.by.noattrs) > 0){
            warning("attributes provided in order.by (",
                    paste(order.by.noattrs, collapse=","),
                    ") are not in attributes and were thus removed." )
        }
        if(length(order.by)==0){
            order.by <- ""
        }else{
            order.by <- paste(order.by, collapse=",")
        }
    }else{
        order.by <- ""
    }
    ## OK, build that query.
    if(order.by!=""){
        orderquery <- paste(" order by", order.by, order.type)
    }else{
        orderquery <- ""
    }
    ## now build the join query that joins all required tables.
    joinquery <- joinQueryOnAttributes(x, attrs=attrs)
    finalquery <- paste0("select distinct ",
                         paste(resultattrs, collapse=","),
                         " from ",
                         joinquery,
                         wherequery,
                         orderquery
                         )
    return(finalquery)
}

.getWhat <- function(x, attrs=listAttributes(x), filter, order.by="",
                     order.type="asc", match.case=FALSE, force=FALSE){
    Q <- .buildQuery(x, attrs=attrs, filter=filter, order.by=order.by,
                     order.type=order.type, match.case=match.case, force=force)
    return(dbGetQuery(dbconn(x), Q))
}


## No need to join tables, as in its present form the database is a single-table database
joinQueryOnAttributes <- function(x, attrs, join, start.table){
    return("mirtarbase")
}

