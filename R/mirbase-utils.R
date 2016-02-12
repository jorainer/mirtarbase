
mirbase.con <- function(){
    return(mirbase_dbconn())
}

mirbase.version <- function(){
    ##    return(dbmeta(mirbase.db:::datacache, "MIRBASESOURCEVERSION"))
    tmp <- mirbase_dbInfo()
    return(tmp[ tmp$name=="MIRBASESOURCEVERSION", "value" ])
}


##*************************************************
##
##  mirna: pre-miRNA information
##  mirna_pre_mature: maps pre-miRNA to mature miRNA on _id (in mirna)
##                    and auto_mature (in mirna_mature)
##  mirna_mature: mature miRNA information
##
##*************************************************
##
##  mirbase.db conversion functions
##
##*************************************************
## maps pre-miRNA names to mature miRNA names (e.g. hsa-mir-16-1 to hsa-miR-16-1-3p and hsa-miR-16-5p)
## x... the pre-miRNA name(s)
## condition... the condition to look for entries
## returns a data.frame with the pre miRNA name and the corresponding mature miRNA name, or NA if not found.
premirna2matmirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2mat(x, condition=condition, ifnotfound=ifnotfound,
                   "mirna_id", "mature_name", return.type=return.type)
    return(Res)
}
premirna2matmirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition,
                   ifnotfound=ifnotfound, pre="mirna_id", mat="mature_acc",
                   return.type=return.type))
}
premirnaAcc2matmirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mirna_acc", mat="mature_acc", return.type=return.type))
}
premirnaAcc2matmirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mirna_acc", mat="mature_name", return.type=return.type))
}
premirnaAcc2premirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mirna_acc", mat="mirna_id", return.type=return.type))
}
premirna2premirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mirna_id", mat="mirna_acc", return.type=return.type))
}
## and the other way round... just switching pre with mat...
matmirna2premirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2mat(x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_name", mat="mirna_id", return.type=return.type)
    return(Res)
}
matmirna2premirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_name", mat="mirna_acc", return.type=return.type))
}
matmirnaAcc2premirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_acc", mat="mirna_acc", return.type=return.type))
}
matmirnaAcc2premirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_acc", mat="mirna_id", return.type=return.type))
}
matmirnaAcc2matmirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_acc", mat="mature_name", return.type=return.type))
}
matmirna2matmirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    return(pre2mat(x=x, condition=condition, ifnotfound=ifnotfound,
                   pre="mature_name", mat="mature_acc", return.type=return.type))
}

##
## CAVE: If we really want to use "like" and ignore.case:
## o The ifnotfound mapping has to be revised.
## o The re-ordering of the results has to be adapted.
pre2mat <- function(x, condition="=", ifnotfound=NA, pre="mirna_id",
                    mat="mature_name", return.type="data.frame", ignore.case=FALSE){
    if(ignore.case){
        nocase <- " collate nocase"
    }else{
        nocase <- ""
    }
    ## well, better to use SQL queries...
    ## straight forward way would be to use the get("hsa-mir-16-1", mirbaseMATURE)
    condition <- match.arg(condition, c("=", "!=", "in", "like", "not like"))
    attrs <- c("mirna_id", "mirna_acc", "mature_name", "mature_acc")
    pre <- match.arg(pre, attrs)
    mat <- match.arg(mat, attrs)
    return.type <- match.arg(return.type, c("data.frame", "list"))
    x.orig <- x
    x <- unique(x)
    if(length(x) > 1){
        if(condition=="=")
            condition <- "in"
        if(condition=="!=")
            condition <- "not in"
        if(!(condition %in% c("in", "not in")))
            stop("condition ", condition, " not allowed! only \"in\" and \"not in\" are valid!")
        x <- paste0("(", paste(sQuote(x), collapse=","), ")")
    }else{
        x <- sQuote(x)
    }
    ## if both pre and mat contain mature I can use a simpler, faster query.
    if(length(grep(c(pre, mat), pattern="^mature")) == 2){
        Q <- paste0("select distinct ", pre, ", ", mat," from mirna_mature where ",
                    pre," ", condition, " ", x)
    }else if(length(grep(c(pre, mat), pattern="^mirna")) == 2){
        Q <- paste0("select distinct ", pre, ", ", mat," from mirna where ",
                    pre," ", condition, " ", x)
    }else{
        Q <- paste0("select distinct ", pre, ", ", mat," from mirna join",
                    " mirna_pre_mature on (mirna._id=mirna_pre_mature._id)",
                    " join mirna_mature on ",
                    "(mirna_pre_mature.auto_mature=mirna_mature.auto_mature)",
                    " where ", pre," ", condition, " ", x, nocase)
    }
    Res <- dbGetQuery(mirbase.con(), Q)
    return(Res)
    ## add NA rows for x if not found:
    notfound <- unique(x.orig[ !(x.orig %in% Res[ , pre ]) ])
    if(length(notfound) > 0){
        tmp <- data.frame(notfound, rep(ifnotfound, length(notfound)))
        colnames(tmp) <- c(pre, mat)
        Res <- rbind(Res, tmp)
    }
    ## want to preserve the input ordering...
    ##idx <- match(Res[ , pre ], x.orig)
    ##Res <- Res[ order(idx), ]
    ##idx <- match(x.orig, Res[ , pre ])
    ##Res <- Res[ idx, ]
    ##if(return.type=="list")
    ##    Res <- split(Res[ , mat ], f=Res[ , pre ])
    ##return(Res)
    ## didn't work... try this:
    if(return.type=="list"){
        Res <- split(Res[ , mat ], Res[ , pre ])
        Res <- Res[ x.orig ]
        return(Res)
    }else{
        Res <- split(Res, Res[ , pre ])
        ## order
        Res <- Res[ x.orig ]
        Res <- do.call(rbind, Res)
        rownames(Res) <- NULL
        return(Res)
    }
}

##************************************************************
##
##   pre-miRNA to mirfam mappings
##
##************************************************************
premirna2mirfam <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mirna_id", "prefam_id", return.type=return.type)
    return(Res)
}
premirnaAcc2mirfam <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mirna_acc", "prefam_id", return.type=return.type)
    return(Res)
}
premirnaAcc2mirfamAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mirna_acc", "prefam_acc", return.type=return.type)
    return(Res)
}
premirna2mirfamAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mirna_id", "prefam_acc", return.type=return.type)
    return(Res)
}
mirfam2premirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_id", fam="mirna_id", return.type=return.type)
    return(Res)
}
mirfam2premirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_id", fam="mirna_acc", return.type=return.type)
    return(Res)
}
mirfamAcc2premirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_acc", fam="mirna_id", return.type=return.type)
    return(Res)
}
mirfamAcc2premirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_acc", fam="mirna_acc", return.type=return.type)
    return(Res)
}
mirfam2mirfamAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_id", fam="prefam_acc", return.type=return.type)
    return(Res)
}
mirfamAcc2mirfam <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- pre2fam(x, condition=condition, ifnotfound=ifnotfound,
                   pre="prefam_acc", fam="prefam_id", return.type=return.type)
    return(Res)
}

## map pre-miRNA ids/acc to miRNA family
pre2fam <- function(x, condition="=", ifnotfound=NA, pre="mirna_id",
                    fam="prefam_id", return.type="data.frame", ignore.case=FALSE){
    if(ignore.case){
        nocase <- " collate nocase"
    }else{
        nocase <- ""
    }
    ## well, better to use SQL queries...
    ## straight forward way would be to use the get("hsa-mir-16-1", mirbaseMATURE)
    condition <- match.arg(condition, c("=", "!=", "in", "like", "not like"))
    attrs <- c("mirna_id", "mirna_acc", "prefam_id", "prefam_acc")
    pre <- match.arg(pre, attrs)
    fam <- match.arg(fam, attrs)
    return.type <- match.arg(return.type, c("data.frame", "list"))
    x.orig <- x
    x <- unique(x)
    if(length(x) > 1){
        if(condition=="=")
            condition <- "in"
        if(condition=="!=")
            condition <- "not in"
        if(!(condition %in% c("in", "not in")))
            stop("condition ", condition, " not allowed! only \"in\" and \"not in\" are valid!")
        x <- paste0("(", paste(sQuote(x), collapse=","), ")")
    }else{
        x <- sQuote(x)
    }
    ## if both pre and mat contain mature I can use a simpler, faster query.
    if(length(grep(c(pre, fam), pattern="^prefam")) == 2){
        Q <- paste0("select distinct ", pre, ", ", fam," from mirna_prefam where ",
                    pre," ", condition, " ", x)
    }else{
        Q <- paste0("select distinct ", pre, ", ", fam," from mirna join",
                    " mirna_2_prefam on (mirna._id=mirna_2_prefam._id) join",
                    " mirna_prefam on",
                    " (mirna_2_prefam.auto_prefam=mirna_prefam.auto_prefam)",
                    " where ", pre," ", condition, " ", x, nocase)
    }
    Res <- dbGetQuery(mirbase.con(), Q)
    ## add NA rows for x if not found:
    notfound <- unique(x.orig[ !(x.orig %in% Res[ , pre ]) ])
    if(length(notfound) > 0){
        tmp <- data.frame(notfound, rep(ifnotfound, length(notfound)))
        colnames(tmp) <- c(pre, fam)
        Res <- rbind(Res, tmp)
    }
    ## want to preserve the input ordering...
    ##idx <- match(Res[ , pre ], x.orig)
    ##Res <- Res[ order(idx), ]
    ##idx <- match(x.orig, Res[ , pre ])
    ##Res <- Res[ idx, ]
    ## didn't work... try this:
    if(return.type=="list"){
        Res <- split(Res[ , fam ], Res[ , pre ])
        Res <- Res[ x.orig ]
        return(Res)
    }else{
        Res <- split(Res, Res[ , pre ])
        ## order
        Res <- Res[ x.orig ]
        Res <- do.call(rbind, Res)
        rownames(Res) <- NULL
        return(Res)
    }
    ##    if(return.type=="list")
    ##    Res <- split(Res[ , fam ], f=Res[ , pre ])
    ##return(Res)
}


##************************************************************
##
##   mature miRNA to mirfam mappings
##
##************************************************************
matmirna2mirfam <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mature_name", "prefam_id", return.type=return.type)
    return(Res)
}
matmirna2mirfamAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mature_name", "prefam_acc", return.type=return.type)
    return(Res)
}
matmirnaAcc2mirfam <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mature_acc", "prefam_id", return.type=return.type)
    return(Res)
}
matmirnaAcc2mirfamAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "mature_acc", "prefam_acc", return.type=return.type)
    return(Res)
}
mirfam2matmirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   mat="prefam_id", fam="mature_name", return.type=return.type)
    return(Res)
}
mirfam2matmirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "prefam_id", "mature_acc", return.type=return.type)
    return(Res)
}
mirfamAcc2matmirna <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "prefam_acc", "mature_name", return.type=return.type)
    return(Res)
}
mirfamAcc2matmirnaAcc <- function(x, condition="=", ifnotfound=NA, return.type="data.frame"){
    Res <- mat2fam(x, condition=condition, ifnotfound=ifnotfound,
                   "prefam_acc", "mature_acc", return.type=return.type)
    return(Res)
}

## map mature miRNA ids/acc to miRNA family
mat2fam <- function(x, condition="=", ifnotfound=NA, mat="mature_name",
                    fam="prefam_id", return.type="data.frame", ignore.case=FALSE){
    if(ignore.case){
        nocase <- " collate nocase"
    }else{
        nocase <- ""
    }
    ## well, better to use SQL queries...
    ## straight forward way would be to use the get("hsa-mir-16-1", mirbaseMATURE)
    condition <- match.arg(condition, c("=", "!=", "in", "like", "not like"))
    attrs <- c("mature_name", "mature_acc", "prefam_id", "prefam_acc")
    pre <- match.arg(mat, attrs)
    fam <- match.arg(fam, attrs)
    return.type <- match.arg(return.type, c("data.frame", "list"))
    x.orig <- x
    x <- unique(x)
    if(length(x) > 1){
        if(condition=="=")
            condition <- "in"
        if(condition=="!=")
            condition <- "not in"
        if(!(condition %in% c("in", "not in")))
            stop("condition ", condition, " not allowed! only \"in\" and \"not in\" are valid!")
        x <- paste0("(", paste(sQuote(x), collapse=","), ")")
    }else{
        x <- sQuote(x)
    }
    Q <- paste0("select distinct ", mat, ", ", fam," from mirna_mature join",
                " mirna_pre_mature on",
                " (mirna_mature.auto_mature=mirna_pre_mature.auto_mature)",
                " join mirna_2_prefam on (mirna_pre_mature._id=mirna_2_prefam._id)",
                " join mirna_prefam on (mirna_2_prefam.auto_prefam=mirna_prefam.auto_prefam)",
                " where ", mat," ", condition, " ", x, nocase)
    Res <- dbGetQuery(mirbase.con(), Q)
    ## add NA rows for x if not found:
    notfound <- unique(x.orig[ !(x.orig %in% Res[ , mat ]) ])
    if(length(notfound) > 0){
        tmp <- data.frame(notfound, rep(ifnotfound, length(notfound)))
        colnames(tmp) <- c(mat, fam)
        Res <- rbind(Res, tmp)
    }
    ## want to preserve the input ordering...
    ##idx <- match(Res[ , mat ], x.orig)
    ##Res <- Res[ order(idx), ]
    ##idx <- match(x.orig, Res[ , mat ])
    ##Res <- Res[ idx, ]
    ##rownames(Res) <- NULL
    ##if(return.type=="list")
    ##    Res <- split(Res[ , fam ], f=Res[ , mat ])
    ##return(Res)
    ## didn't work... try this:
    if(return.type=="list"){
        Res <- split(Res[ , fam ], Res[ , mat ])
        Res <- Res[ x.orig ]
        return(Res)
    }else{
        Res <- split(Res, Res[ , mat ])
        ## order
        Res <- Res[ x.orig ]
        Res <- do.call(rbind, Res)
        rownames(Res) <- NULL
        return(Res)
    }
}


