## some basic tests for methods of the database.

test_EntrezidFilter <- function(){
    EFilt <- EntrezFilter(123, condition="contains")
    EFilt
    ## default attribute of the filter
    checkEquals(column(EFilt), "entrez")
    checkEquals(where(EFilt), "like '%123%'")
    ## attribute of the filter for the mirtarbase package:
    checkEquals(column(EFilt, mirtarbase), "target_gene_entrez_gene_id")
    ## and the where query:
    checkEquals(where(EFilt, mirtarbase),
                "target_gene_entrez_gene_id like '%123%'")

    ## EntrezidFilter with more than one value:
    EFilt <- EntrezFilter(c(123, 3435))
    EFilt
    checkEquals(where(EFilt, mirtarbase),
                "target_gene_entrez_gene_id in ('123','3435')")
    ## What if we want to set >
    checkException(EntrezFilter("sddsf", condition = ">"))
}

test_ExperimentFilter <- function(){
    ## create a simple experiment filter; this would return all entries that
    ## have "qpcr" in the experiments attribute; like is case insensitive.
    EFilt <- ExperimentFilter("qPCR", condition = "contains")
    checkEquals(column(EFilt, mirtarbase), "experiments")
    checkEquals(where(EFilt, mirtarbase), "experiments like '%qPCR%'")
    checkEquals(where(EFilt), "like '%qPCR%'")
}


test_GenenameFilter <- function(){
    ## create a gene name filter
    GFilt <- GenenameFilter("Bcl2l11", condition = "!=")
    GFilt
    ## the attribute name in the ensembldb package
    checkEquals(column(GFilt), "gene_name")
    ## and in the mirtarbase package
    checkEquals(column(GFilt, mirtarbase), "target_gene")
    checkEquals(where(GFilt, mirtarbase), "target_gene !='Bcl2l11'")
}

test_MatmirnaFilter <- function(){
    ## mature miRNA filter
    MFilt <- MatMirnaFilter(c("hsa-miR-16-5p", "hsa-miR-16-3p"))
    MFilt
    ## and in the mirtarbase package
    checkEquals(column(MFilt, mirtarbase), "mirna")
    checkEquals(where(MFilt, mirtarbase),
                "mirna in ('hsa-miR-16-5p','hsa-miR-16-3p')")
}

test_PublicationFilter <- function(){
    ## create a simple publication filter (on Pubmed IDs)
    PFilt <- PublicationFilter(123)
    PFilt
    ## the attribute in the mirtarbase package
    checkEquals(column(PFilt, mirtarbase), "references_pmid")
    checkEquals(where(PFilt, mirtarbase), "references_pmid ='123'")
}

test_SpeciesFilter <- function(){
    ## create a species filter on genes
    SFilt <- SpeciesFilter("Homo sapiens", feature="gene")
    SFilt
    ## the attribute for the gene
    checkEquals(column(SFilt, mirtarbase), "species_target_gene")
    checkEquals(where(SFilt, mirtarbase), "species_target_gene ='Homo sapiens'")

    ## the same for miRNA
    SFilt <- SpeciesFilter("Homo sapiens", feature="mirna")
    ## the attribute for miRNA
    checkEquals(column(SFilt, mirtarbase), "species_mirna")
    checkEquals(where(SFilt, mirtarbase), "species_mirna ='Homo sapiens'")

    ## what if we used a non-standard species?
    SFilt <- SpeciesFilter("Homo hobbeniensis", feature="gene")
    ## this throws a warning
    suppressWarnings(
        checkEquals(where(SFilt, mirtarbase),
                    "species_target_gene ='Homo hobbeniensis'")
    )

    ## where has to be the same for "Homo sapiens" and "Homo_sapiens",
    ## if we're using a mirtarbase.
    A <- SpeciesFilter("Homo sapiens")
    B <- SpeciesFilter("Homo_sapiens")
    checkTrue(where(A) != where(B))
    checkTrue(where(A, mirtarbase) == where(B, mirtarbase))
}

test_SupportTypeFilter <- function(){
    ## create a filter on the mirtarbase support types
    SFilt <- SupportTypeFilter("unknown")
    SFilt
    ## the attribute that will be queried
    checkEquals(column(SFilt, mirtarbase), "support_type")
    ## throws a warning since the support type unknown is not present in the database
    suppressWarnings(
        checkEquals(where(SFilt, mirtarbase), "support_type ='unknown'")
    )
    SFilt <- SupportTypeFilter("Functional MTI")
    checkEquals(where(SFilt, mirtarbase), "support_type ='Functional MTI'")
}

test_MirtarbaseidFilter <- function(){

    mtf <- MirtarbaseIdFilter("unknown")
    ## the attribute that will be queried
    checkEquals(column(mtf, mirtarbase), "mirtarbase_id")
    checkEquals(where(mtf, mirtarbase), "mirtarbase_id ='unknown'")
}


test_PremirnaFilter <- function(){
    ## pre miRNA filter.
    ## note that the PremirnaFilter methods for MirtarbaseDb internally map the
    ## pre-miRNA name to mature miRNA names.
    PF <- PreMirnaFilter(c("bla", "hsa-mir-15b"))
    checkEquals(where(PF), "in ('bla','hsa-mir-15b')")
    ## doesn't make much sense...
    checkEquals(column(PF, mirtarbase), "mirna")
    ## pre-miRNAs get *translated* automagically into mature miRNA names.
    checkEquals(where(PF, mirtarbase),
                "mirna in ('hsa-miR-15b-5p','hsa-miR-15b-3p')")

    ## what if nothing was found...
    PF <- PreMirnaFilter("baf")
    checkEquals(where(PF, mirtarbase), NULL)
    ## returns nothing

    PF <- PreMirnaFilter(c("bla", "hsa-mir-15b"), condition="!=")
    checkEquals(where(PF, mirtarbase),
                "mirna not in ('hsa-miR-15b-5p','hsa-miR-15b-3p')")

    ##
    PF <- PreMirnaFilter("hsa-mir-15b")
    checkEquals(where(PF, mirtarbase),
                "mirna in ('hsa-miR-15b-5p','hsa-miR-15b-3p')")
}


test_PremirnaidFilter <- function(){
    PF <- PreMirnaIdFilter("MI0000070")
    checkEquals(where(PF, mirtarbase),
                "mirna in ('hsa-miR-16-5p','hsa-miR-16-1-3p')")
}

test_MatmirnaidFilter <- function(){
    MF <- MatMirnaIdFilter("MIMAT0004489", "!=")
    checkEquals(where(MF, mirtarbase), "mirna !='hsa-miR-16-1-3p'")
    MF <- MatMirnaIdFilter(c("MIMAT0004489", "dfd"), "==")
    checkEquals(where(MF, mirtarbase), "mirna ='hsa-miR-16-1-3p'")
}


test_MirfamFilter <- function(){
    MF <- MirfamFilter("mir-15")
    checkTrue(length(where(MF, mirtarbase)) > 0)
}

test_MirfamidFilter <- function(){
    MF <- MirfamIdFilter("MIPF0000006")
    checkTrue(length(where(MF, mirtarbase)) > 0)
}




