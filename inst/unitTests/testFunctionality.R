test_mtis <- function(){
    ## get all MTIs for gene BCL2 and restrict to human.
    BCL2 <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2"),
                                         SpeciesFilter("Homo sapiens", feature="gene")))
    checkEquals(unique(gene(BCL2)), "BCL2")
    checkTrue(any(matmirna(BCL2) == "hsa-miR-16-5p"))
    ## for how many miRNAs there is evidence that they target this gene?
    length(BCL2)
    ## just looking at the first
    BCL2[[1]]

    ## all of em
    BCL2

    ## what evidences are there for the interaction?
    table(unlist(lapply(BCL2, supportedBy)))

    ## which pre-miRNAs are these?
    pres <- lapply(BCL2, premirna)
    pres2 <- premirna(BCL2)
    checkEquals(pres, pres2)
    BCL2[ which(is.na(pres)) ]
    ## well, some of the mature miRNAs can not be mapped to pre-miRNAs, because
    ## the mirbase version between the mirbase.db and the miRTarbase does not
    ## match.

    ## Check if it's working with a single Filter
    BCL2 <- mtis(mirtarbase, filter=GenenameFilter("BCL2"))

    ## Check what happens if we use return.type="data.frame"
    Test <- mtis(mirtarbase, filter=GenenameFilter("BCL2"), return.type="data.frame")
    checkTrue(nrow(Test) != length(unique(Test$target_gene)))
    Test <- mtis(mirtarbase, columns=c("target_gene"), filter=GenenameFilter("BCL2"),
                 return.type="data.frame")
    checkTrue(nrow(Test) == length(unique(Test$target_gene)))
    Test <- mtis(mirtarbase, columns=c("mirna", "target_gene"),
                 filter=list(GenenameFilter("BCL2"), SpeciesFilter("Homo sapiens")),
                 return.type="data.frame")
    checkEquals(nrow(Test), length(unique(Test$mirna)))

    ## Is there any "enrichment" for a specific miRNA family?
    tmp <- mirfam(BCL2)
    tmp2 <- unlist(lapply(BCL2, mirfam))
    checkEquals(tmp, tmp2)
    checkEquals(names(sort(table(tmp), decreasing=TRUE))[1], "mir-15")

    ## now we perform the same query as above, but ask for a data.frame as
    ## return type.
    BCL2.df <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2"), SpeciesFilter("Homo sapiens", feature="gene")), return.type="data.frame")

    ## the query is much faster, but:
    nrow(BCL2.df)
    length(BCL2)
    ## we get almost twice as many rows: one row in the data.frame corresponds to one
    ## publication in which a miRNA target gene interaction was detected, while each
    ## element in the list corresponds to one MTI class with the list of Report objects
    ## (i.e. publications) being stored in the report slot.

    ## Just testing what happens if we query for a non-existing pre-miRNA
    checkException(mtis(mirtarbase, filter=list(PremirnaFilter("hsa-mir-181z"),
                                                SpeciesFilter("Homo sapiens"))))

    ## Seems to be a problem with the mirbase.db version mismatch.
    ## REPORT THAT IN THE HELP AND IN THE VIGNETTE
    Test <- mtis(mirtarbase, filter=list(PremirnaFilter("hsa-mir-181d"),
                                         SpeciesFilter("Homo sapiens")))
    premirna2matmirna("hsa-mir-181d")
    Test <- mtis(mirtarbase, filter=list(MatmirnaFilter("hsa-miR-181d%", condition="like"),
                                         SpeciesFilter("Homo sapiens")))
    unique(matmirna(Test))
    ## Thus, apparently, the mature miRNA has been renamed from hsa-miR-181d to
    ## "hsa-miR-181d-5p".

    ## get all MTIs for a specific miRNA
    Test <- mtis(mirtarbase, filter=list(PremirnaFilter(c("hsa-mir-223")),
                                         SpeciesFilter("Homo sapiens")))
    length(Test)
    ## get the gene names along with the number of supporting evidences
    do.call(rbind, lapply(Test, function(x){ return(c(gene=gene(x),
                                                      report_count=reportCount(x))) }))

    ## get all MTIs for a miRNA family.
    MTIs <- mtis(mirtarbase, filter=list(MirfamFilter("mir-15"),
                                         SpeciesFilter("Homo sapiens")))
    length(MTIs)
    ## next we might want to ask if there are genes targeted by more than one miRNA
    ## of this family.
    head(sort(table(gene(MTIs)), decreasing=TRUE))
}

test_premirna2matmirna <- function(){
    ## detach("package:mirtarbase", unload=TRUE)
    ## library(mirtarbase)
    tmp <- premirna2matmirna("hsa-miR-181%", condition="like")
    ## For like we want to get many rows.
    checkTrue(nrow(tmp) > 1)
    tmp2 <- premirna2matmirna("hsa-mir-181%", condition="like")
    ## And it shouldn't matter wheter we're case sensitive or not.
    checkEquals(tmp, tmp2)

    ## We don't want to have the pattern in the results too.
    checkTrue(!any(tmp$mirna_id == "hsa-miR-181%"))

    ## If we don't find anything, we expect to get NA
    tmp <- premirna2matmirna("asdfdsf", condition="like")
    checkEquals(as.character(tmp[1, 1]), "asdfdsf")

    ## LLLL go on here.
    ## What else to test:
    ## 1) a pre-miRNA encoding 2 mature miRNAs
    tmp <- premirna2matmirna("hsa-miR-16-2")  ## dont' expect anything here.
    checkEquals(tmp[1, 2], NA)
    tmp <- premirna2matmirna("hsa-mir-16-2")
    checkEquals(tmp[, 2], c("hsa-miR-16-5p", "hsa-miR-16-2-3p"))
    ## 2) two pre-miRNAs encoding 3 mature miRNAs
    tmp <- premirna2matmirna(c("hsa-mir-16-1", "hsa-mir-16-2"))
    checkEquals(tmp[, 2], c("hsa-miR-16-5p", "hsa-miR-16-1-3p", "hsa-miR-16-5p",
                            "hsa-miR-16-2-3p"))

    ## Assume we're using case insensitive queries here!
    premirna2matmirna("hsa-mir-181a-2")
    premirna2matmirna("hsa-miR-181a-2")

}

test_pre2fam <- function(){
    ## detach("package:mirtarbase", unload=TRUE)
    ## library(mirtarbase)
    tmp <- premirna2mirfam("hsa-mir-15%", condition="like")
    tmp2 <- premirna2mirfam("hsa-miR-15%", condition="like")
    checkEquals(tmp, tmp2)
    checkTrue(!any(tmp$mirna_id == "hsa-mir-15%"))

    tmp <- premirna2mirfam("%-mir-15b", condition="like")
    tmp2 <- premirna2mirfam("hsa-mir-15b")
    checkEquals(unique(tmp[, "prefam_id"]), unique(tmp2[, "prefam_id"]))
    ##
    ## Check if ordering is preserved.
    tmp <- premirna2mirfam(c("hsa-mir-16-2", "hsa-mir-15b", "hsa-mir-16-2"))
    checkEquals(tmp[, 1], c("hsa-mir-16-2", "hsa-mir-15b", "hsa-mir-16-2"))

    ## Other way round.
    tmp <- mirfam2premirna("mir-15")
    checkTrue(any(tmp[, 2] == "hsa-mir-15b"))
}

test_mat2fam <- function(){
    ## detach("package:mirtarbase", unload=TRUE)
    ## library(mirtarbase)
    tmp <- matmirna2mirfam("hsa-miR-15%", condition="like")
    tmp2 <- matmirna2mirfam("hsa-mir-15%", condition="like")
    checkEquals(tmp, tmp2)
    checkTrue(!any(tmp[, 1] == "hsa-miR-15%"))

    tmp <- matmirna2mirfam("%-miR-15b-5p", condition="like")
    tmp2 <- matmirna2mirfam("hsa-miR-15b-5p")
    checkEquals(unique(tmp[, 2]), unique(tmp2[, 2]))
    ##
    ## Check if ordering is preserved.
    tmp <- matmirna2mirfam(c("hsa-miR-16-3p", "hsa-miR-15b-5p", "hsa-miR-16-3p"))
    checkEquals(tmp[, 1], c("hsa-miR-16-3p", "hsa-miR-15b-5p", "hsa-miR-16-3p"))

    ## Other way round.
    tmp <- mirfam2matmirna("mir-15")
    checkTrue(any(tmp[, 2] == "hsa-miR-15b-5p"))
}



##*****************
##
## test some stuff with MTIList
##*****************
test_mtilist <- function(x){
    detach("package:mirtarbase", unload=TRUE)
    library(mirtarbase)

    ## get all MTIs for a miRNA family.
    MTIs <- mtis(mirtarbase, filter=list(PremirnaFilter("hsa-mir-223"),
                                         SpeciesFilter("Homo sapiens")))
    checkTrue(is(MTIs, "MTIList"))

    fullDf <- (as.data.frame(MTIs, collapse.reports=NULL))
    collDf <- as.data.frame(MTIs, collapse.reports=";")
    checkTrue(nrow(fullDf) != nrow(collDf))

    head(entrezid(MTIs))

    head(experiments(MTIs))

    head(gene(MTIs))

    head(geneSpecies(MTIs))

    head(id(MTIs))

    head(matmirna(MTIs))

    head(mirnaSpecies(MTIs))

    head(pmid(MTIs))

    head(supportedBy(MTIs))

    head(premirna(MTIs))

    head(reportCount(MTIs))

    head(mirfam(MTIs))


    ## get an MTI from the database we're looking here for an interaction
    ## between the gene BCL2L11 and miRNAs from the mir-17 family.
    BCL2L11 <- mtis(mirtarbase, filter=list(GenenameFilter("BCL2L11"), MirfamFilter("mir-17")))
    BCL2L11 <- BCL2L11[[1]]
    BCL2L11

    ## what's the mature miRNA?
    matmirna(BCL2L11)
    ## what's the precursor of this miRNA
    premirna(BCL2L11)

    ## the gene name and the Entrezgene ID?
    gene(BCL2L11)
    entrezid(BCL2L11)

    ## on how many reports (publications) does this base?
    reportCount(BCL2L11)

    ## what evidence/support types?
    supportedBy(BCL2L11)

    ## what experiments:
    experiments(BCL2L11)

    ## list the reports:
    reports(BCL2L11)
}

##******************************************************
##
##    the mtisBy
##
##******************************************************
## for comparison...
test_mtisBy <- function(){
    Filters <- list(PremirnaFilter(c("hsa-mir-15b", "hsa-mir-16-2", "hsa-mir-17")),
                    SpeciesFilter("Homo sapiens"),
                    GenenameFilter(c("BCL2", "BCL2L11", "MCL1")))
    Test <- mtis(mirtarbase, filter=Filters, return.type="MTI")
    length(Test)
    Test

    ## MTIs by gene
    TestBy <- mtisBy(mirtarbase, by="gene", filter=Filters)
    TestBy

    checkEquals(length(TestBy), 3)
    ## MTIs by support type
    TestBy <- mtisBy(mirtarbase, by="support_type", filter=Filters)
    TestBy

    ## MTIs by mirfam
    TestBy <- mtisBy(mirtarbase, by="mirfam", filter=Filters)
    TestBy

    ## MTIs by premirna
    TestBy <- mtisBy(mirtarbase, by="premirna", filter=Filters)
    TestBy


    ## MTIs by publication
    TestBy <- mtisBy(mirtarbase, by="pmid", filter=Filters)
    TestBy

    ## works. what if we get a little more: can we run by premirna on more?
    TestBy <- mtisBy(mirtarbase, by="premirna", filter=list(MirfamFilter(c("mir-17", "mir-15")),
                                                            SpeciesFilter("Homo sapiens")))
    TestBy

}


