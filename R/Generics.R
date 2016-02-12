##******************************************************
##
##  Generic methods
##
##******************************************************
if(!isGeneric("as.data.frame")){
    setGeneric("as.data.frame", function(x, ...)
        standardGeneric("as.data.frame"))
}
if(!isGeneric("attribute")){
    setGeneric("attribute", function(object, db, ...)
        standardGeneric("attribute"))
}
## if(!isGeneric("connection")){
##     setGeneric("connection", function(x, ...)
##         standardGeneric("connection"))
## }
if(!isGeneric("entrezid"))
    setGeneric("entrezid" , function(object, ...)
               standardGeneric("entrezid"))
if(!isGeneric("experiments"))
    setGeneric("experiments", function(object, ...)
        standardGeneric("experiments"))
if(!isGeneric("gene"))
    setGeneric("gene" , function(object, ...)
        standardGeneric("gene"))
if(!isGeneric("geneSpecies"))
    setGeneric("geneSpecies" , function(object, ...)
               standardGeneric("geneSpecies"))
if(!isGeneric("id"))
    setGeneric("id" , function(object, ...)
               standardGeneric("id"))
if(!isGeneric("listAttributes")){
    setGeneric("listAttributes", function(x, ...)
        standardGeneric("listAttributes"))
}
if(!isGeneric("listExperiments")){
    setGeneric("listExperiments", function(x, ...)
        standardGeneric("listExperiments"))
}
if(!isGeneric("listPmids")){
    setGeneric("listPmids", function(x, ...)
        standardGeneric("listPmids"))
}
if(!isGeneric("listSpecies")){
    setGeneric("listSpecies", function(x, ...)
        standardGeneric("listSpecies"))
}
if(!isGeneric("listSupportTypes")){
    setGeneric("listSupportTypes", function(x, ...)
        standardGeneric("listSupportTypes"))
}
if(!isGeneric("matmirna"))
    setGeneric("matmirna", function(object, ...)
        standardGeneric("matmirna"))
if(!isGeneric("matmirnaId"))
    setGeneric("matmirnaId" , function(object)
               standardGeneric("matmirnaId"))
## if(!isGeneric("matmirnaSequence"))
##     setGeneric("matmirnaSequence" , function(object, ...)
##                standardGeneric("matmirnaSequence"))
if(!isGeneric("mirfam"))
    setGeneric("mirfam" , function(object, ...)
               standardGeneric("mirfam"))
if(!isGeneric("mirnaSpecies"))
    setGeneric("mirnaSpecies", function(object, ...)
        standardGeneric("mirnaSpecies"))
if(!isGeneric("mtis"))
    setGeneric("mtis", function(x, ...)
        standardGeneric("mtis"))
if(!isGeneric("mtisBy"))
    setGeneric("mtisBy", function(x, by, ...)
        standardGeneric("mtisBy"))
if(!isGeneric("pmid"))
    setGeneric("pmid", function(object, ...)
               standardGeneric("pmid"))
if(!isGeneric("premirna"))
    setGeneric("premirna" , function(object, ...)
               standardGeneric("premirna"))
if(!isGeneric("premirnaId"))
    setGeneric("premirnaId" , function(object)
               standardGeneric("premirnaId"))
setGeneric("requireTable", function(x, db, ...)
    standardGeneric("requireTable"))
if(!isGeneric("reports<-"))
    setGeneric("reports<-", function(object, value)
        standardGeneric("reports<-"))
if(!isGeneric("reportCount"))
    setGeneric("reportCount" , function(object, ...)
               standardGeneric("reportCount"))
if(!isGeneric("reports"))
    setGeneric("reports", function(x, ...)
        standardGeneric("reports"))
if(!isGeneric("shortShow"))
    setGeneric("shortShow", function(object, ...)
        standardGeneric("shortShow"))
if(!isGeneric("supportedBy"))
    setGeneric("supportedBy", function(object, ...)
        standardGeneric("supportedBy"))
if(!isGeneric("listTables"))
    setGeneric("listTables", function(x, ...)
        standardGeneric("listTables"))
if(!isGeneric("version"))
    setGeneric("version", function(object, ...)
        standardGeneric("version"))

