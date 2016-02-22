## some basic tests for methods of the database.

test_EntrezidFilter <- function(){
    ##***********************
    ##
    ## EntrezidFilter
    ##***********************
    ##
    EFilt <- EntrezidFilter(123, condition="like")
    EFilt
    ## default attribute of the filter
    checkEquals(column(EFilt), "entrezid")
    ## attribute of the filter for the mirtarbase package:
    checkEquals(column(EFilt, mirtarbase), "target_gene_entrez_gene_id")
    ## and the where query:
    where(EFilt, mirtarbase)

    ## EntrezidFilter with more than one value:
    EFilt <- EntrezidFilter(c(123, 3435))
    EFilt
    where(EFilt, mirtarbase)
    ## What if we want to set >
    checkException(EntrezidFilter("sddsf", condition = ">"))
}

test_ExperimentFilter <- function(){
    ##***********************
    ##
    ## ExperimentFilter
    ##***********************
    ## create a simple experiment filter; this would return all entries that have "qpcr"
    ## in the experiments attribute; like is case insensitive.
    EFilt <- ExperimentFilter("%qPCR%", condition="like")
    checkEquals(column(EFilt, mirtarbase), "experiments")
    where(EFilt, mirtarbase)
    ## like with more than two values; not that like does not support multiple values,
    ## thus it is reverted to "in"
    checkException(EFilt <- ExperimentFilter(c("qPCR", "microarray"), condition="like"))

    ## a good starting point is to list all experiments from the database
    ## and select the pattern from there
    listExperiments(mirtarbase)
}


test_GenenameFilter <- function(){
    ##***********************
    ##
    ## GenenameFilter
    ##***********************
    ## create a gene name filter
    GFilt <- GenenameFilter("Bcl2l11", condition="!=")
    GFilt
    ## the attribute name in the ensembldb package
    checkEquals(column(GFilt), "gene_name")
    ## and in the mirtarbase package
    checkEquals(column(GFilt, mirtarbase), "target_gene")
    where(GFilt, mirtarbase)
}

test_MatmirnaFilter <- function(){
    ##***********************
    ##
    ## MatmirnaFilter
    ##***********************
    ## mature miRNA filter
    MFilt <- MatmirnaFilter(c("hsa-miR-16-5p", "hsa-miR-16-3p"))
    MFilt
    ## and in the mirtarbase package
    checkEquals(column(MFilt, mirtarbase), "mirna")
    where(MFilt, mirtarbase)
}

test_PublicationFilter <- function(){
    ##***********************
    ##
    ## PublicationFilter
    ##***********************
    ## create a simple publication filter (on Pubmed IDs)
    PFilt <- PublicationFilter(123)
    PFilt
    ## the attribute in the mirtarbase package
    checkEquals(column(PFilt, mirtarbase), "references_pmid")
    where(PFilt, mirtarbase)
}

test_SpeciesFilter <- function(){
    ##***********************
    ##
    ## SpeciesFilter
    ##***********************
    ## create a species filter on genes
    SFilt <- SpeciesFilter("Homo sapiens", feature="gene")
    SFilt
    ## the attribute for the gene
    checkEquals(column(SFilt, mirtarbase), "species_target_gene")
    where(SFilt, mirtarbase)

    ## the same for miRNA
    SFilt <- SpeciesFilter("Homo sapiens", feature="mirna")
    ## the attribute for miRNA
    checkEquals(column(SFilt, mirtarbase), "species_mirna")
    where(SFilt, mirtarbase)

    ## what if we used a non-standard species?
    SFilt <- SpeciesFilter("Homo hobbeniensis", feature="gene")
    ## this throws a warning
    suppressWarnings(
        where(SFilt, mirtarbase)
    )

    ## where has to be the same for "Homo sapiens" and "Homo_sapiens",
    ## if we're using a mirtarbase.
    A <- SpeciesFilter("Homo sapiens")
    B <- SpeciesFilter("Homo_sapiens")
    checkTrue(where(A) != where(B))
    checkTrue(where(A, mirtarbase) == where(B, mirtarbase))

    ## obviously only values present in the database make sense:
    listSpecies(mirtarbase, "gene")
    listSpecies(mirtarbase, "mirna")
}

test_SupportTypeFilter <- function(){

    ##***********************
    ##
    ## SupportTypeFilter
    ##***********************
    ## create a filter on the mirtarbase support types
    SFilt <- SupportTypeFilter("unknown")
    SFilt
    ## the attribute that will be queried
    checkEquals(column(SFilt, mirtarbase), "support_type")
    ## throws a warning since the support type unknown is not present in the database
    suppressWarnings(
        where(SFilt, mirtarbase)
    )

    ## better to select one of the correct support types:
    listSupportTypes(mirtarbase)
}

test_MirtarbaseidFilter <- function(){

    mtf <- MirtarbaseidFilter("unknown")
    ## the attribute that will be queried
    checkEquals(column(mtf, mirtarbase), "mirtarbase_id")
    where(mtf, mirtarbase)

}


test_PremirnaFilter <- function(){
    ##***********************
    ##
    ## PremirnaFilter
    ##***********************
    ## pre miRNA filter.
    ## note that the PremirnaFilter methods for MirtarbaseDb internally map the pre-miRNA name
    ## to mature miRNA names.
    PF <- PremirnaFilter(c("bla", "hsa-mir-15b"))
    where(PF)
    ## doesn't make much sense...
    column(PF, mirtarbase)  ## that's the database column in which we're going to search
    ## the where internally maps the pre-miRNA ID to the corresponding mature miRNA ID(s)
    where(PF, mirtarbase)

    ## what if nothing was found...
    PF <- PremirnaFilter("baf")
    where(PF, mirtarbase)
    ## returns nothing

    PF <- PremirnaFilter(c("bla", "hsa-mir-15b"), condition="!=")
    where(PF, mirtarbase)

    ##
    PF <- PremirnaFilter("hsa-mir-15b")
    where(PF , mirtarbase)
}


test_PremirnaidFilter <- function(){
    ##***********************
    ##
    ## PremirnaidFilter
    ##***********************
    PF <- PremirnaidFilter("MI0000070")
    where(PF, mirtarbase)
}

test_MatmirnaidFilter <- function(){
    ##***********************
    ##
    ## MatmirnaidFilter
    ##***********************
    MF <- MatmirnaidFilter("MIMAT0004489", "!=")
    where(MF, mirtarbase)


    MF <- MatmirnaidFilter(c("MIMAT0004489", "dfd"), "in")
    where(MF, mirtarbase)
}


test_MirfamFilter <- function(){
    ##***********************
    ##
    ## MirfamFilter
    ##***********************
    MF <- MirfamFilter("mir-15")
    where(MF, mirtarbase)
}

test_MirfamidFilter <- function(){
    ##***********************
    ##
    ## MirfamidFilter
    ##***********************
    MF <- MirfamidFilter("MIPF0000006")
    where(MF, mirtarbase)
}




