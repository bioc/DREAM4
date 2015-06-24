library(RUnit)
library(SummarizedExperiment)
#-------------------------------------------------------------------------------
runTests <- function()
{
    stopifnot(file.exists("insilico_size10_1"))
    stopifnot(file.exists("insilico_size100_1"))
    
    dir10 <- file.path("insilico_size10_1")
    dir100 <- file.path("insilico_size100_1")

    test_buildDream4SummarizedExperimentSize_10(dir10)
    test_buildDream4SummarizedExperimentSize_100(dir100)

} # runTests
#-------------------------------------------------------------------------------
buildDream4SummarizedExperiment <- function(sourceDirectory)
{
    file.ts <- file.path(sourceDirectory, 'timeseries.tsv')
    file.kd <- file.path(sourceDirectory, 'knockdowns.tsv')
    file.ko <- file.path(sourceDirectory, 'knockouts.tsv')
    file.mf <- file.path(sourceDirectory, 'multifactorial.tsv')
    file.wt <- file.path(sourceDirectory, 'wildtype.tsv')
    file.dk <- file.path(sourceDirectory, 'dualknockouts.tsv')
    file.gs <- file.path(sourceDirectory, 'goldStandard.tsv')


    stopifnot(file.exists(file.ts))
    tbl.ts <- read.table(file.ts, sep='\t', header=TRUE)

    stopifnot(file.exists(file.kd))
    tbl.kd <- read.table(file.kd, sep='\t', header=TRUE)

    stopifnot(file.exists(file.ko))
    tbl.ko <- read.table(file.ko, sep='\t', header=TRUE)

    stopifnot(file.exists(file.mf))
    tbl.mf <- read.table(file.mf, sep='\t', header=TRUE)

    stopifnot(file.exists(file.wt))
    tbl.wt <- read.table(file.wt, sep='\t', header=TRUE)

    stopifnot(file.exists(file.dk))
    tbl.dk <- read.table(file.dk, sep='\t', header=T, as.is<-T)

    stopifnot(file.exists(file.gs))
    tbl.gs = read.table(file.gs, sep='\t', header=T, as.is=T)

    tbl <- t(tbl.wt)
    colnames(tbl) [1] <- 'wt'

     # The files *timeseries.tsv contain time courses showing how the network
     # responds to a perturbation and how it relaxes upon removal of the
     # perturbation. For networks of size 10 we provide 5 different time
     # series, for networks of size 100 we provide 10 time series. Each time
     # series has 21 time points. The initial condition always corresponds to a
     # steady-state measurement of the wild-type. At t=0, a perturbation is
     # applied to the network as described below. The first half of the time
     # series (until t=500) shows the response of the network to the
     # perturbation. At t=500, the perturbation is removed (the wild-type
     # network is restored). The second half of the time series (until t=1000)
     # shows how the gene expression levels go back from the perturbed to the
     # wild-type state.  In contrast to the multifactorial perturbations
     # described in the previous section, which affect all the genes
     # simultaneously, the perturbations applied here only affect about a third
     # of all genes, but basal activation of these genes can be strongly
     # increased or decreased. For example, these experiments could correspond
     # to physical or chemical perturbations applied to the cells, which would
     # cause (via regulatory mechanisms not explicitly modeled here) some genes
     # to have an increased or decreased basal activation. The genes that are
     # directly targeted by the perturbation may then cause a change in the
     # expression level of their downstream target genes.
   
     # the timeseries table has Time in the first column 21 time points, from 0
     # 1000 by increments of 50 there is further complexity: DREAM specifies
     # t0: wild type t50-450: the perturbation (ie, oxygen) is applied t500
     # G1-Gmax make up the remainder of the columns since we are building an
     # expression matrix where the rows are genes, and the various conditions
     # are columns, after transposing, we must pick out the

    x <- t(tbl.ts [, 2:ncol(tbl.ts)])
    new.colnames <- as.character(tbl.ts$Time)
    rep.factor <- 10 # number of timesteps
  
    prefixes <- c('perturbation.1',  rep('perturbation.1.applied',  10),
                  rep('perturbation.1.removed',  10),
                  'perturbation.2',  rep('perturbation.2.applied',  10),
                  rep('perturbation.2.removed',  10), 
                  'perturbation.3',  rep('perturbation.3.applied',  10),
                  rep('perturbation.3.removed',  10), 
                  'perturbation.4',  rep('perturbation.4.applied',  10),
                  rep('perturbation.4.removed',  10), 
                  'perturbation.5',  rep('perturbation.5.applied',  10),
                  rep('perturbation.5.removed',  10),
                  'perturbation.6',  rep('perturbation.6.applied',  10),
                  rep('perturbation.6.removed',  10),
                  'perturbation.7',  rep('perturbation.7.applied',  10),
                  rep('perturbation.7.removed',  10), 
                  'perturbation.8',  rep('perturbation.8.applied',  10),
                  rep('perturbation.8.removed',  10), 
                  'perturbation.9',  rep('perturbation.9.applied',  10),
                  rep('perturbation.9.removed',  10), 
                  'perturbation.10', rep('perturbation.10.applied', 10),
                  rep('perturbation.10.removed', 10))
      
    if(nrow(tbl) == 10)
        prefixes <- prefixes [1:105]
    else if(nrow(tbl) == 100)
        prefixes <- prefixes [1:210]

    new.colnames <- paste(prefixes, new.colnames, sep='.t')
    colnames(x) <- new.colnames
    tbl <- cbind(tbl, x)         

    x <- t(tbl.ko)
    colnames(x) <- paste(rownames(x), '.ko', sep='')
    tbl <- cbind(tbl, x)

    
    x <- t(tbl.kd)
    colnames(x) <- paste(rownames(x), '.kd', sep='')
    tbl <- cbind(tbl, x)


      # The files *multifactorial.tsv contain steady-state levels of variations
      # of the network, which are obtained by applying multifactorial
      # perturbations to the original network. Each line gives the steady state
      # of a different perturbation experiment, i.e., of a different variation
      # of the network. One may think of each experiment as a gene expression
      # profile from a different patient, for example. We simulate
      # multifactorial perturbations by slightly increasing or decreasing the
      # basal activation of all genes of the network simultaneously by
      # different random amounts.

    x <- t(tbl.mf)
    column.count <- ncol(x)
    colnames(x) <- paste('MF.', 1:column.count, sep='')
    tbl <- cbind(tbl, x)

    sumexp <- SummarizedExperiment(SimpleList(simulated=tbl))

      # the goldstandard matrix 
    gs <- matrix(rep(0,(column.count * column.count)), nrow=column.count)

    rownames(gs) <- rownames(tbl)
    colnames(gs) <- rownames(tbl)
    x <- tbl.gs [which(tbl.gs [,3] == 1),]
    for(r in 1:nrow(x))
      gs [x$G1[r], x$G2[r]] <- 1

    metadata(sumexp) <- list(goldStandardAdjacencyMatrix=gs)

    metadata(sumexp) [['doubleKnockoutGenePairs']] <- tbl.dk

    invisible(sumexp)

} # build.dream4.rdata.file 
#-------------------------------------------------------------------------------
test_buildDream4SummarizedExperimentSize_10 <- function(sourceDir)
{
    print('--- test.build.dream4.summarizedExperiment.size.10')

    se <- buildDream4SummarizedExperiment(sourceDirectory=sourceDir)

    #checkEquals(sort(slotNames(se)), c("assays", "colData", "metadata",
    #                                   "rowRanges"))

    mtx.expr <- assays(se)[[1]]
    checkEquals(dim(mtx.expr), c(10, 136))  

    column.data <- colData(se)   # column meta data for possibly multiple
                                 # assay matrices are all shared
    row.data <- rowRanges(se)      # row meta data for possibly multiple assay
                                 # matrices are all shared
    experiment.data <- metadata(se)

    checkEquals(nrow(column.data), ncol(mtx.expr))
    checkEquals(names(row.data), rownames(mtx.expr))

    checkEquals(names(experiment.data),
                c("goldStandardAdjacencyMatrix", "doubleKnockoutGenePairs"))

} # test_buildDream4SummarizedExperimentSize_10
#-------------------------------------------------------------------------------
test_buildDream4SummarizedExperimentSize_100 <- function(sourceDir)
{
    print('--- test.build.dream4.summarizedExperiment.sizd.100')
    se <- buildDream4SummarizedExperiment(sourceDirectory=sourceDir)
    #checkEquals(sort(slotNames(se)), c("assays", "colData", "metadata",
    #                                   "rowRanges"))
    mtx.expr <- assays(se)[[1]]
    checkEquals(dim(mtx.expr), c(100, 511))
      # check these columns.  there should be 21 for each perturbation,
      # and 10 perturbations total for the size 100 network:  210 total
    checkEquals(length(grep('perturbation', colnames(mtx.expr))), 210)

    column.data <- colData(se)   # column meta data for possibly multiple
                                 # assay matrices are all shared
    row.data <- rowRanges(se)      # row meta data for possibly multiple assay
                                 # matrices are all shared
    experiment.data <- metadata(se)
  
    checkEquals(nrow(column.data), ncol(mtx.expr))
    checkEquals(names(row.data), rownames(mtx.expr))

    checkEquals(names(experiment.data),
                c("goldStandardAdjacencyMatrix", "doubleKnockoutGenePairs"))

} # test_buildDream4SummarizedExperimentSize_100
#-------------------------------------------------------------------------------
run <- function(dataDir, destDir)
{
    challenges <- c(dream4_010_01='insilico_size10_1',
                     dream4_010_02='insilico_size10_2',
                     dream4_010_03='insilico_size10_3',
                     dream4_010_04='insilico_size10_4',
                     dream4_010_05='insilico_size10_5',
                     dream4_100_01='insilico_size100_1',
                     dream4_100_02='insilico_size100_2',
                     dream4_100_03='insilico_size100_3',
                     dream4_100_04='insilico_size100_4',
                     dream4_100_05='insilico_size100_5')

    for(name in names(challenges)) {
        cmd <-
          sprintf("%s = buildDream4SummarizedExperiment(sourceDirectory='%s')",
                  name, challenges [[name]])
        print(cmd)
        eval(parse(text<-cmd))
        filename <- sprintf('%s/%s.RData', destDir, name)
        cmd <- sprintf('save(%s, file="%s")', name, filename)
        print(cmd)
        eval(parse(text=cmd))
        }
  
} # run
#-------------------------------------------------------------------------------
