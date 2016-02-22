####============================================================
##  Test functions for AnnotationDbi methods
##
####------------------------------------------------------------
## detach("package:mirtarbase", unload=TRUE)
## library(mirtarbase)

test_columns <- function(){
    cols <- columns(mirtarbase)
    checkTrue(any(cols == "ENTREZID"))
    checkTrue(!any(cols == "not there"))

    ## Check mappings
    mapped <- mirtarbase:::.mapCols2Mirtarbase(c("ENTREZID", "MATMIRNA",
                                                 "PMID"))
    checkEquals(mapped, c("target_gene_entrez_gene_id", "mirna", "references_pmid"))
    ## Check exception
    checkException(mirtarbase:::.mapCols2Mirtarbase(c("ENTREZID", "PREMIRNANAME")))

    ## Reverse mapping
    mapped <- mirtarbase:::.mapMirtarbase2Cols(c("species_target_gene", "mirtarbase_id",
                                                 "support_type"))
    checkEquals(mapped, c("GENESPECIES", "MIRTARBASEID", "SUPPORTTYPE"))
    ## Check exception
    checkException(mirtarbase:::.mapMirtarbase2Cols(c("mirtarbase_id", "don't exist")))
}

test_filter4cols <- function(){
    checkException(mirtarbase:::filterForKeytype("PREMIRNA"))

    ## Check some filters.
    filts <- mirtarbase:::filterForKeytype("PMID")
    checkTrue(is(filts, "PublicationFilter"))
    filts <- mirtarbase:::filterForKeytype("GENESPECIES")
    checkTrue(is(filts, "SpeciesFilter"))
    checkTrue(filts@feature == "gene")
    filts <- mirtarbase:::filterForKeytype("MIRNASPECIES")
    checkTrue(is(filts, "SpeciesFilter"))
    checkTrue(filts@feature == "mirna")
    ## Other filters
    checkTrue(is(mirtarbase:::filterForKeytype("SUPPORTTYPE"), "SupportTypeFilter"))
    checkTrue(is(mirtarbase:::filterForKeytype("EXPERIMENT"), "ExperimentFilter"))
    checkTrue(is(mirtarbase:::filterForKeytype("MIRTARBASEID"), "MirtarbaseidFilter"))
    ## Filters defined in ensembldb
    checkTrue(is(mirtarbase:::filterForKeytype("SYMBOL"), "GenenameFilter"))
    checkTrue(is(mirtarbase:::filterForKeytype("ENTREZID"), "EntrezidFilter"))
    ## Filters defined in mirhostgenes
    checkTrue(is(mirtarbase:::filterForKeytype("MATMIRNA"), "MatmirnaFilter"))

}

test_keys <- function(){
    ks <- keys(mirtarbase)
    checkTrue(length(grep(ks[1], pattern="^MIRT")) > 0)

    ## Check different keytypes.
    ks <- keys(mirtarbase, keytype="ENTREZID")
    checkTrue(is.numeric(ks))

    ks <- keys(mirtarbase, keytype="GENESPECIES")
    checkEquals(length(unique(ks)), length(ks))
    ks <- keys(mirtarbase, keytype="MIRNASPECIES")
    checkEquals(length(unique(ks)), length(ks))

    ks <- keys(mirtarbase, keytype="MATMIRNA")
    ks <- keys(mirtarbase, keytype="PMID")
    checkEquals(length(unique(ks)), length(ks))

    ks <- keys(mirtarbase, keytype="SUPPORTTYPE")
    checkEquals(ks, listSupportTypes(mirtarbase))
    ks <- keys(mirtarbase, keytype="SYMBOL")
    ## Check exception.
    checkException(keys(mirtarbase, keytype="BLA"))

    ## And what with experiments???
    ks <- keys(mirtarbase, keytype="EXPERIMENT")
}

test_select <- function(){
    ## Select gene symbol.
    Test <- select(mirtarbase, keys="BCL2", keytype="SYMBOL")
    checkTrue(all(colnames(Test) %in% columns(mirtarbase)))

    Test2 <- select(mirtarbase, keys=GenenameFilter("BCL2"))
    checkEquals(Test, Test2)

    ## Error check:
    checkException(select(mirtarbase, keys="BCL2"))
    checkException(select(mirtarbase, keys="BCL2", keytype="bo"))

    ## Select mature miRNA.
    Test <- select(mirtarbase, keys="hsa-miR-15b-5p", keytype="MATMIRNA")
    checkTrue(all(Test$MATMIRNA == "hsa-miR-15b-5p"))
    ## Combine
    Test2 <- select(mirtarbase, keys=list(MatmirnaFilter("hsa-miR-15b-5p"),
                                          SupportTypeFilter("Functional MTI")))
    Test <- Test[Test$SUPPORTTYPE == "Functional MTI", ]
    Test <- Test[order(Test$MIRTARBASEID),]
    Test2 <- Test2[order(Test2$MIRTARBASEID), ]
    rownames(Test) <- NULL
    rownames(Test2) <- NULL
    checkEquals(Test2, Test)

    ## Select support type
    Test <- select(mirtarbase, keys="Functional MTI", keytype="SUPPORTTYPE")
    Test2 <- select(mirtarbase, keys=SupportTypeFilter("Functional MTI"))
    checkEquals(Test, Test2)
    Test <- select(mirtarbase, keys=list(SupportTypeFilter("Functional MTI"),
                                         SpeciesFilter("Homo sapiens", feature="mirna")))
    checkTrue(all(Test$SUPPORTTYPE == "Functional MTI"))
    checkTrue(all(Test$MIRNASPECIES == "Homo sapiens"))

}



