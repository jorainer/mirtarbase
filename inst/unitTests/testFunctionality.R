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
}

test_premirna2matmirna <- function(){
    premirna2matmirna("hsa-miR-181%", condition="like")
    premirna2matmirna("hsa-mir-181%", condition="like")

    ## Assume we're using case insensitive queries here!
    premirna2matmirna("hsa-mir-181a-2")
    premirna2matmirna("hsa-miR-181a-2")

    premirna2matmirna("hsa-mir-16-2")

}

## is there any "enrichment" for a specific miRNA family?
sort(table(unlist(lapply(BCL2, mirfam))), decreasing=TRUE)

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

##\dontrun{
## just testing what happens if we query for a pre-miRNA
Test <- mtis(mirtarbase, filter=list(PremirnaFilter("hsa-mir-181z"), SpeciesFilter("Homo sapiens")))
## throws an error.
##}

## seems to be a problem with the mirbase.db version mismatch
Test <- mtis(mirtarbase, filter=list(PremirnaFilter("hsa-mir-181d"), SpeciesFilter("Homo sapiens")))

## get all MTIs for a specific miRNA
Test <- mtis(mirtarbase, filter=list(PremirnaFilter(c("hsa-mir-223")), SpeciesFilter("Homo sapiens")))
length(Test)
## get the gene names along with the number of supporting evidences
do.call(rbind, lapply(Test, function(x){ return(c(gene=gene(x), report_count=reportCount(x))) }))

## get all MTIs for a miRNA family.
MTIs <- mtis(mirtarbase, filter=list(MirfamFilter("mir-15"), SpeciesFilter("Homo sapiens")))
length(MTIs)
## next we might want to ask if there are genes targeted by more than one miRNA
## of this family.
head(sort(table(unlist(lapply(MTIs, gene))), decreasing=TRUE))

##*****************
##
## test some stuff with MTIList
##*****************
MTIL <- MTIList(MTIs)

MTIL

head(as.data.frame(MTIL))

head(entrezid(MTIL))

head(experiments(MTIL))

head(gene(MTIL))

head(geneSpecies(MTIL))

head(id(MTIL))

head(matmirna(MTIL))

head(mirnaSpecies(MTIL))

head(pmid(MTIL))

head(supportedBy(MTIL))

head(premirna(MTIL))

head(reportCount(MTIL))

head(mirfam(MTIL))


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


##******************************************************
##
##    the mtisBy
##
##******************************************************
## for comparison...
Filters <- list(PremirnaFilter(c("hsa-mir-15b", "hsa-mir-16-2", "hsa-mir-17")),
                SpeciesFilter("Homo sapiens"),
                GenenameFilter(c("BCL2", "BCL2L11", "MCL1")))
Test <- mtis(mirtarbase, filter=Filters, return.type="MTI")
length(Test)
Test

## MTIs by gene
TestBy <- mtisBy(mirtarbase, by="gene", filter=Filters)
TestBy

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
TestBy <- mtisBy(mirtarbase, by="pmid", filter=list(SpeciesFilter("Homo sapiens")))
head(TestBy)

## cross check that.
Res <- dbGetQuery(connection(mirtarbase), "select * from mirtarbase where species_target_gene='Homo sapiens'")
## should be OK.

## works. what if we get a little more: can we run by premirna on more?
TestBy <- mtisBy(mirtarbase, by="premirna", filter=list(MirfamFilter(c("mir-17", "mir-15"))))
length(TestBy)
TestBy

