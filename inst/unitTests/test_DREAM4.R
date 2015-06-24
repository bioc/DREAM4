library(DREAM4)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests = function()
{
  test_loadSmallDataSet()
  test_loadLargeDataSet()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_loadSmallDataSet = function()
{
    print("--- loadSmallDataSet")
    data(dream4_010_01)  
    #checkEquals(sort(slotNames(dream4_010_01)), c("assays", "colData", "metadata", "rowRanges"))
    checkEquals(nrow(assays(dream4_010_01)[[1]]), 10)

         # do some sanity checks on the expression matrix
    mtx <- assays(dream4_010_01)[[1]]
    checkEquals(dim(mtx), c(10, 136))
    geneNames <- paste("G", 1:10, sep="")
    checkEquals(rownames(mtx), geneNames)
    checkEquals(head(colnames(mtx)),
                c("wt", "perturbation.1.t0", "perturbation.1.applied.t50",
                  "perturbation.1.applied.t100", "perturbation.1.applied.t150",
                  "perturbation.1.applied.t200"))

        # and on the accompanying data
    checkEquals(names(metadata(dream4_010_01)),
                c("goldStandardAdjacencyMatrix", "doubleKnockoutGenePairs"))

    mtx.goldStandard <- metadata(dream4_010_01)[[1]]
    checkEquals(rownames(mtx.goldStandard), geneNames)
    checkEquals(colnames(mtx.goldStandard), geneNames)
    
    dbl.ko <- metadata(dream4_010_01)[[2]]
    checkEquals(dim(dbl.ko), c(5,2))
    checkEquals(colnames(dbl.ko), c ("G_i", "G_j"))

        # metadata on the expression matrix: an empty GRanges for each row
    checkEquals(names(rowRanges(dream4_010_01)), geneNames)
    checkEquals(length(rowRanges(dream4_010_01)[["G1"]]), 0)
        # no metadata for the columns
    checkEquals(names(colData(dream4_010_01)), character(0))


} # test_loadSmallDataSet
#------------------------------------------------------------------------------------------------------------------------
test_loadLargeDataSet = function()
{
    print("--- loadLargeDataSet")
    data(dream4_100_05)  
    #checkEquals(sort(slotNames(dream4_100_05)), c("assays", "colData", "metadata", "rowRanges"))
    checkEquals(nrow(assays(dream4_100_05)[[1]]), 100)

        # do some sanity checks on the expression matrix
    mtx <- assays(dream4_100_05)[[1]]
    checkEquals(dim(mtx), c(100, 511))
    geneNames <- paste("G", 1:100, sep="")
    checkEquals(rownames(mtx), geneNames)
    
    checkEquals(head(colnames(mtx)),
                c("wt", "perturbation.1.t0", "perturbation.1.applied.t50",
                  "perturbation.1.applied.t100", "perturbation.1.applied.t150",
                  "perturbation.1.applied.t200"))
    checkEquals(length(grep("perturbation", colnames(mtx))), 210)

       # and now on the accompanying data
    checkEquals(names(metadata(dream4_100_05)),
                c("goldStandardAdjacencyMatrix", "doubleKnockoutGenePairs"))

    mtx.goldStandard <- metadata(dream4_100_05)[[1]]
    checkEquals(rownames(mtx.goldStandard), geneNames)
    checkEquals(colnames(mtx.goldStandard), geneNames)
    
    dbl.ko <- metadata(dream4_100_05)[[2]]
    checkEquals(dim(dbl.ko), c(20,2))
    checkEquals(colnames(dbl.ko), c ("G_i", "G_j"))

        # metadata on the expression matrix: an empty GRanges for each row
    checkEquals(names(rowRanges(dream4_100_05)), geneNames)
    checkEquals(length(rowRanges(dream4_100_05)[["G1"]]), 0)
        # no metadata for the columns
    checkEquals(names(colData(dream4_100_05)), character(0))


} # test_loadLargeDataSet
#------------------------------------------------------------------------------------------------------------------------
