## some basic tests for methods of the database.
detachem <- function( x ){
    NS <- loadedNamespaces()
    if( any( NS==x ) ){
        pkgn <- paste0( "package:", x )
        detach( pkgn, unload=TRUE, character.only=TRUE )
    }
}
Pkgs <- c( "mirtarbase" )
tmp <- sapply( Pkgs, detachem )
tmp <- sapply( Pkgs, library, character.only=TRUE )


Q <- mirtarbase:::.buildQuery( mirtarbase )
Q

## can we use regexp???
## no way! they would have to be registered using the sqlite3_create_function which doesn't seem to be accessible through RSQLite.
## example; http://www.blackdogfoundry.com/blog/supporting-regular-expressions-in-sqlite/
dbGetQuery( connection( mirtarbase ), "select * from mirtarbase limit 3" )
dbGetQuery( connection( mirtarbase ), "select * from mirtarbase where target_gene regexp '^DUS'" )



SF <- SpeciesFilter( "Mus musculus" )
SF
Q <- mirtarbase:::.buildQuery( mirtarbase, filter=list( SF ) )
Q
tmp <- dbGetQuery( connection( mirtarbase ), Q )
head( tmp )
nrow( tmp )


## do something more specific: get all for hsa-miR-16-5p.
F <- MatmirnaFilter( "hsa-miR-16-5p" )
Res <- mirtarbase:::.getWhat( mirtarbase, filter=list( F ) )
nrow( Res )

## well ... way too many...
F <- GenenameFilter( "BCL2L11" )
Res <- mirtarbase:::.getWhat( mirtarbase, filter=list( F ) )
nrow( Res )
Res
## as we see, the query is performed case insensitive! If we don't want that
## we have to say match.case=TRUE
Res <- mirtarbase:::.getWhat( mirtarbase, filter=list( F ), match.case=TRUE )
nrow( Res )
Res
## or restrict to human genes only would also do it...
Res <- mirtarbase:::.getWhat( mirtarbase, filter=list( F, SpeciesFilter( "Homo sapiens", feature="gene" ) ), order.by="mirna" )
nrow( Res )
Res
## we get now however duplicated miRNA-target gene pairs, one for each publication
## in which the MTI was validated


## check the performance of mclapply vs normal apply...
## 10 rows
Res <- dbGetQuery( connection( mirtarbase ), "select * from mirtarbase order by mirtarbase_id limit 10" )
system.time(
    bli <- mirtarbase:::data.frame2mti( Res )
    )
## 0.008
system.time(
    bla <- mirtarbase:::data.frame2mtiNreport( Res )
    )
## 0.016
## check if it makes sense:
Res[ Res$mirtarbase_id=="MIRT000006", ]
bla[["MIRT000006"]]
## Perfectly fine!!!


## 50 rows:
Res <- dbGetQuery( connection( mirtarbase ), "select * from mirtarbase order by mirtarbase_id limit 50" )
system.time(
    bli <- mirtarbase:::data.frame2mti( Res )
    )
## 0.024
system.time(
    bla <- mirtarbase:::data.frame2mtiNreport( Res )
    )
## 0.077


## 100 rows:
Res <- dbGetQuery( connection( mirtarbase ), "select * from mirtarbase order by mirtarbase_id limit 100" )
system.time(
    bli <- mirtarbase:::data.frame2mti( Res )
    )
## 0.039
system.time(
    bla <- mirtarbase:::data.frame2mtiNreport( Res )
    )
## 0.151


## 400 rows:
Res <- dbGetQuery( connection( mirtarbase ), "select * from mirtarbase order by mirtarbase_id limit 400" )
system.time(
    bli <- mirtarbase:::data.frame2mti( Res )
    )
## 0.150
system.time(
    bla <- mirtarbase:::data.frame2mtiNreport( Res )
    )
## 0.526


## 1000 rows:
Res <- dbGetQuery( connection( mirtarbase ), "select * from mirtarbase order by mirtarbase_id limit 1000" )
system.time(
    bli <- mirtarbase:::data.frame2mti( Res )
    )
## 0.395
system.time(
    bla <- mirtarbase:::data.frame2mtiNreport( Res )
    )
## 1.369

Test <- bli[[1]]



##*********************************************************
##
##  Testing conversion functions
##
##*********************************************************
## pre-miRNA to mature miRNA
Res <- premirna2matmirna( c( "bla", "hsa-mir-16-1", "hsa-mir-15b", "hsa-mir-16-2" ), return.type="list" )
Res
Res <- premirna2matmirna( c( "bla", "hsa-mir-16-1", "hsa-mir-15b", "hsa-mir-16-2", "hsa-mir-16-1" ) )
Res
Res <- premirna2matmirnaAcc( c( "hsa-mir-15b" ) )
Res
Res <- premirnaAcc2matmirna( c( "MI0000070" ) )
Res
Res <- premirnaAcc2matmirnaAcc( c( "MI0000070" ) )
Res
Res <- premirna2premirnaAcc( "hsa-mir-15b" )
Res
Res <- premirnaAcc2premirna( "MI0000070" )
Res

## mature miRNA to pre-miRNA
Res <- matmirna2premirna( c( "hsa-miR-16-5p", "hsa-miR-16-3p", "hsa-miR-16-1-3p" ) )
Res
Res <- matmirna2premirnaAcc( c( "hsa-miR-16-5p", "hsa-miR-16-3p" ) )
Res
Res <- matmirnaAcc2premirna( c( "MIMAT0000069", "MIMAT0004489" ) )
Res
Res <- matmirnaAcc2premirnaAcc( "MIMAT0004489" )
Res
Res <- matmirna2matmirnaAcc( "hsa-miR-16-5p" )
Res
Res <- matmirnaAcc2matmirna( "MIMAT0004489" )
Res

## pre-miRNA to mirfam and vice versa
Res <- premirna2mirfam( c( "hsa-mir-16-1", "hsa-mir-15b", "hsa-mir-16-2", "hsa-mir-15b" ) )
Res
Res <- premirnaAcc2mirfam( "MI0000070" )
Res
Res <- mirfam2premirna( "mir-15" )
head( Res )
Res <- mirfam2premirnaAcc( "mir-15" )
head( Res )
Res <- mirfam2mirfamAcc( "mir-15" )
Res
Res <- mirfamAcc2mirfam( "MIPF0000006" )
Res
Res <- mirfamAcc2premirna( "MIPF0000006" )
head( Res )


## mature miRNAs to mirfam and vice versa
Res <- matmirna2mirfam( "hsa-miR-15b-5p" )
Res
Res <- matmirna2mirfamAcc( c( "bla", "hsa-miR-15b-5p" ) )
Res
Res <- mirfam2matmirna( "mir-15" )   ## weird to get only a handful of mature miRNAs...
head( Res )
Res <- mirfamAcc2matmirna( "MIPF0000006" )
head( Res )
Res <- mirfam2matmirnaAcc( "mir-15" )
head( Res )

