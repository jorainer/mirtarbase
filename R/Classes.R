## Definition of all classes.

##*************************************************
##
##  MirtarbaseDb
##
##  the main class representing the database connection.
##
##*************************************************
setClass("MirtarbaseDb",
         representation(con="DBIConnection",
                        tables="list",
                        mirtarbase_version="character",
                        mirtarbase_date="character",
                        species_target_gene="character",
                        species_mirna="character",
                        support_type="character"
                        ),
         prototype=list(con=NULL,
                        tables=list(),
                        mirtarbase_version="",
                        mirtarbase_date="",
                        species_target_gene="",
                        species_mirna="",
                        support_type=""
                        )
         )




##***********************************************************************
##
##      Main classes representing the data in the database.
##
##***********************************************************************
## A report (i.e. publication) that reported the MTI.
setClass("MTIReport",
         representation(pmid="numeric",
                        experiments="character",
                        support_type="character"
                        ),
         prototype=list(pmid=0,
                        experiments="",
                        support_type=""
                        )
         )
MTIReport <- function(pmid=0, support_type="", experiments=""){
    if(is.na(pmid)){
        return(NA)
    }
    rep <- new("MTIReport",
               pmid=pmid,
               support_type=support_type,
               experiments=experiments)
    return(rep)
}

## the main MTI class representing a miRNA target gene interaction.
setClass("MTI",
         representation(id = "character",
                        mature_mirna="character",
                        species_mirna="character",
                        query="character",
                        target_gene="character",
                        target_gene_entrezid="numeric",
                        species_target_gene="character",
                        report="list"),
         prototype=list(id="",
                        mature_mirna="",
                        species_mirna="",
                        query="",
                        target_gene="",
                        target_gene_entrezid=0,
                        species_target_gene="",
                        report=list(new("MTIReport"))
                        )
         )
MTI <- function(id="", mature_mirna="", species_mirna="", query="", target_gene="",
                target_gene_entrezid=0, species_target_gene="", report=list()){
    ## a MTI has to have an ID!
    if(is.na(id)){
        return(NA)
    }
    mti <- new("MTI",
               id=id,
               mature_mirna=mature_mirna,
               species_mirna=species_mirna,
               query=query,
               target_gene=target_gene,
               target_gene_entrezid=target_gene_entrezid,
               species_target_gene=species_target_gene,
               report=report
               )
    return(mti)
}

## the MTIList class
setClass("MTIList",
         contains="SimpleList",
         representation(),
         prototype(elementType="MTI", listData=list())
         )
## constructor.
MTIList <- function(...){
    args <- list(...)
    if(length(args)==1L && is.list(args[[1L]]))
        args <- args[[1L]]


    new_SimpleList_from_list("MTIList", args)
}

## that's from the S4Vectors package:
new_SimpleList_from_list <- function(Class, x, ..., mcols)
{
    if (!extends(Class, "SimpleList"))
        stop("class ", Class, " must extend SimpleList")
    if (!is.list(x))
        stop("'x' must be a list")
    if (is.array(x)) { # drop any unwanted dimensions
        tmp_names <- names(x)
        dim(x) <- NULL # clears the names
        names(x) <- tmp_names
    }
    class(x) <- "list"
    ans_elementType <- elementType(new(Class))
    if (!all(sapply(x, function(xi) extends(class(xi), ans_elementType))))
        stop("all elements in 'x' must be ", ans_elementType, " objects")
    if (missing(mcols))
        return(new2(Class, listData=x, ..., check=FALSE))
    new2(Class, listData=x, ..., elementMetadata=mcols, check=FALSE)
}


##***********************************************************************
##
##     Filter classes
##
##     all Filters form ensembldb are supported.
##     We're defining there in addition:
##
##     SpeciesFilter
##     PublicationFilter
##     ExperimentFilter
##     SupportTypeFilter
##
##***********************************************************************
## filter for species
setClass("ExperimentFilter", contains="CharacterFilter",
         prototype = list(
             condition = "==",
             value="",
             field = "experiment"
         )
         )
ExperimentFilter <- function(value, condition = "=="){
    return(new("ExperimentFilter", condition = condition,
               value=as.character(value)))
}

## filter for species
setClass("PublicationFilter", contains="CharacterFilter",
         prototype=list(
             condition = "==",
             value = "",
             field = "pmid"
         )
         )
PublicationFilter <- function(value, condition = "=="){
    return(new("PublicationFilter", condition=condition,
               value=as.character(value)))
}

## filter for species
setClass("SpeciesFilter", contains="CharacterFilter",
         representation(feature="character"),
         prototype=list(
             feature = "gene",
             condition = "==",
             value="",
             field = "species"
         )
         )
SpeciesFilter <- function(value, condition = "==", feature = "gene"){
    feature <- match.arg(feature, c("gene", "mirna"))
    return(new("SpeciesFilter", condition=condition,
               value=as.character(value), feature=feature))
}

## filter for species
setClass("SupportTypeFilter", contains="CharacterFilter",
         prototype=list(
             condition = "==",
             value = "",
             field = "support_type"
         )
         )
SupportTypeFilter <- function(value, condition = "=="){
    return(new("SupportTypeFilter", condition=condition,
               value=as.character(value)))
}

## filter for species
setClass("MirtarbaseIdFilter", contains="CharacterFilter",
         prototype=list(
             condition = "==",
             value = "",
             field = "mirtarbase_id"
         )
         )
MirtarbaseIdFilter <- function(value, condition = "=="){
    return(new("MirtarbaseIdFilter", condition=condition,
               value=as.character(value)))
}

