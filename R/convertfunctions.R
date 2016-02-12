### functions that convert data.frames to objects and objects to data.frames

######
##
## Report
##
## converts a data.frame into a list of Report classes (length is nrow(data.frame))
## calls unique on the column names required to build a Report object to avoid having duplicate entries in the result.
data.frame2report <- function(x, colname.pmid="references_pmid",
                              colname.support="support_type",
                              colname.experiments="experiments", split="//", do.unique=TRUE){
    ## check for require columns
    Need <- c(colname.pmid, colname.support, colname.experiments)
    if(class(x)!="data.frame"){
        stop("Can only take data.frame as input!")
    }
    if(sum(colnames(x) %in% Need)!=length(Need)){
        stop(paste0("One or more required columns (", paste(Need, collapse=","), ") not found!"))
    }
    if(do.unique){
        x <- unique(x[ , Need, drop=FALSE ])
    }
    ## so, now I'm happy.
    return(apply(x, MARGIN=1, row2report, name.pmid=colname.pmid, name.support=colname.support, name.experiments=colname.experiments, split=split))
}

## create a (single row) data.frame for a Report object
report2data.frame <- function(x, colname.pmid="references_pmid",
                              colname.support="support_type",
                              colname.experiments="experiments",
                              split="//", stringsAsFactors=FALSE, row.names=NULL){
    DF <- data.frame(a=x@pmid, b=paste(x@experiments, collapse=split),
                     c=x@support_type, stringsAsFactors=stringsAsFactors, row.names=row.names)
    colnames(DF) <- c(colname.pmid, colname.experiments, colname.support)
    return(DF)
}


######
##
## MTI
##
## converts a data.frame into a list of MTI classes (length is nrow(data.frame))
## Note: we are NOT loading the report slot here!!!
data.frame2mti <- function(x, colname.id="mirtarbase_id",
                           colname.mirna="mirna", colname.mirna.species="species_mirna",
                           colname.target="target_gene",
                           colname.target.entrezid="target_gene_entrez_gene_id",
                           colname.target.species="species_target_gene", do.unique=TRUE){
    Need <- c(colname.id, colname.mirna, colname.mirna.species, colname.target,
              colname.target.entrezid, colname.target.species)
    if(class(x)!="data.frame"){
        stop("Can only take data.frame as input!")
    }
    if(sum(colnames(x) %in% Need)!=length(Need)){
        stop(paste0("One or more required columns (", paste(Need, collapse=","), ") not found!"))
    }
    if(do.unique){
        x <- unique(x[ , Need, drop=FALSE ])
    }
    return(apply(x, MARGIN=1, row2mti, name.id=colname.id,
                 name.mirna=colname.mirna, name.mirna.species=colname.mirna.species,
                 name.target=colname.target, name.target.entrezid=colname.target.entrezid,
                 name.target.species=colname.target.species))
}

## create a data.frame for a MTI object, the number of rows depending on the number of reports.
## if collapse.reports not NULL, the reports will be collapsed (separated by collapse.reports).
mti2data.frame <- function(x, colname.id="mirtarbase_id", colname.mirna="mirna",
                           colname.mirna.species="species_mirna",
                           colname.target="target_gene",
                           colname.target.entrezid="target_gene_entrez_gene_id",
                           colname.target.species="species_target_gene",
                           collapse.reports=NULL, stringsAsFactors=FALSE, row.names=NULL){
    ## do first the Report stuff...
    reps <- reports(x)
    reps.df <- do.call("rbind", lapply(reps, report2data.frame,
                                       stringsAsFactors=stringsAsFactors))
    ## OK, now I have a data.frame with nrow=length(report).
    ## collapsing the reports if wished...
    if(!is.null(collapse.reports)){
        reps.df <- data.frame(sapply(reps.df, paste,
                                     collapse=collapse.reports,
                                     simplify=FALSE), stringsAsFactors = stringsAsFactors)
    }
    ## defining the data.frame for the MTI.
    mti.df <- data.frame(a=x@id, b=x@mature_mirna, c=x@species_mirna,
                         d=x@target_gene, e=x@target_gene_entrezid,
                         f=x@species_target_gene, stringsAsFactors=stringsAsFactors,
                         row.names=row.names)
    colnames(mti.df) <- c(colname.id, colname.mirna, colname.mirna.species,
                          colname.target, colname.target.entrezid, colname.target.species)
    ## replicate rows in mti.df to match nrow reps.df
    mti.df <- mti.df[ rep(1, each=nrow(reps.df)), ]
    mti.df <- cbind(mti.df, report_count=reportCount(x), reps.df,
                    stringsAsFactors=stringsAsFactors)
    return(mti.df)
}


## split the data.frame by mirtarbase id.
## (mc)lapply over the list and produce both in one go, MTI from the first row, report from the remaining.
data.frame2mtiNreport <- function(x, colname.pmid="references_pmid",
                                  colname.support="support_type",
                                  colname.experiments="experiments", split="//",
                                  colname.id="mirtarbase_id", colname.mirna="mirna",
                                  colname.mirna.species="species_mirna",
                                  colname.target="target_gene",
                                  colname.target.entrezid="target_gene_entrez_gene_id",
                                  colname.target.species="species_target_gene"){
    x.split <- split(x, f=x[ , colname.id ])
    FUN <- lapply
    ## we can expect a speedup if there are more than 100 MTIs to process...
    if(length(x.split) > 100)
        FUN <- mclapply
    Res <- FUN(x.split, function(z){
        ## build the Report.
        Reps <- apply(z, MARGIN=1,
                      row2report,
                      name.pmid=colname.pmid,
                      name.support=colname.support,
                      name.experiments=colname.experiments,
                      split=split
                     )
        ## build the MTI
        mti <- MTI(id=z[ 1, colname.id ],
                   mature_mirna=z[ 1, colname.mirna ],
                   species_mirna=z[ 1, colname.mirna.species ],
                   target_gene=z[ 1, colname.target ],
                   target_gene_entrezid=z[ 1, colname.target.entrezid ],
                   species_target_gene=z[ 1, colname.target.species],
                   report=Reps
                  )
        return(mti)
    })

    return(MTIList(Res))
}


##############
### internal functions!
row2report <- function(x, name.pmid="references_pmid", name.support="support_type",
                       name.experiments="experiments", split="//"){
    ## here we do not check again for available data etc.
    Rp <- MTIReport(pmid=as.numeric(x[ name.pmid ]),
                 support_type=x[ name.support ],
                 experiments=unique(unlist(strsplit(x[ name.experiments ], split=split)))
                )
    return(Rp)
}
row2mti <- function(x, name.id="mirtarbase_id", name.mirna="mirna",
                    name.mirna.species="species_mirna", name.target="target_gene",
                    name.target.entrezid="target_gene_entrez_gene_id",
                    name.target.species="species_target_gene"){
    ## here we do not check again for available data etc.
    mti <- MTI(
        id=x[ name.id ],
        mature_mirna=x[ name.mirna ],
        species_mirna=x[ name.mirna.species ],
        target_gene=x[ name.target ],
        target_gene_entrezid=as.numeric(x[ name.target.entrezid ]),
        species_target_gene=x[ name.target.species ]
   )
    return(mti)
}

