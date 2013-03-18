### R code from vignette source 'DREAM4.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: init
###################################################
printf <- function(...)print(noquote(sprintf(...)))
options(width=60)


###################################################
### code chunk number 2: load-01-01
###################################################
library(DREAM4)
data(dream4_010_01)
show(dream4_010_01)


###################################################
### code chunk number 3: extractFullMatrix
###################################################
mtx.all <- assays (dream4_010_01)[[1]]
dim(mtx.all)
mtx.all[1:10, c(1,2,126)]
set.seed(37)
print(colnames(mtx.all)[sort(sample(1:ncol(mtx.all), 16))])


###################################################
### code chunk number 4: retrieveExpressionData
###################################################
mtx.goldStandard <- exptData (dream4_010_01)[[1]]


###################################################
### code chunk number 5: colnames
###################################################
grep("perturbation.1.", colnames(mtx.all), v=T, fixed=TRUE)


###################################################
### code chunk number 6: extractAndTranspose
###################################################
ts1.columns <- grep("perturbation.1.", colnames(mtx.all), fixed=TRUE)
mtx.ts1 <- t(mtx.all[, ts1.columns])
print(mtx.ts1[1:3, 1:3])


###################################################
### code chunk number 7: extractGoldStandard
###################################################
print(names(exptData(dream4_010_01)))
mtx.goldStandard <- exptData(dream4_010_01)[[1]]


###################################################
### code chunk number 8: transformGS
###################################################
idx <- which(mtx.goldStandard == 1)
idx.m1 <- idx -1
rows <- idx.m1 %% nrow (mtx.goldStandard) + 1
cols <- idx.m1 %/% nrow (mtx.goldStandard) + 1
tbl.goldStandard <- data.frame(Regulator=rownames(mtx.goldStandard)[rows],
                             Target=colnames(mtx.goldStandard)[cols],
                             Source=rep('goldStandard', length(rows)),
                             stringsAsFactors=FALSE)


###################################################
### code chunk number 9: infer0
###################################################
library(networkBMA)
tbl.inferred <- networkBMA(mtx.ts1, nTimePoints=nrow(mtx.ts1), self=FALSE)
tbl.inferred <- tbl.inferred[order(tbl.inferred$PostProb, decreasing=TRUE),]


###################################################
### code chunk number 10: aupr0
###################################################
tbl.contingency <- contabs.netwBMA(tbl.inferred, tbl.goldStandard[,-3])
#precision <- with(tbl.contingency, TP/(TP+FP))
#recall <- with(tbl.contingency, TP/(TP+FN))
pr <- scores(tbl.contingency, what = c("precision","recall"))
colors <- c ("blue", "darkred")
plot(pr$recall, pr$precision, type='b',
     xlab='RECALL', ylab='PRECISION', 
     col=colors[1],
     xlim=c(0,1), ylim=c(0,1))


###################################################
### code chunk number 11: printaupr
###################################################
print(prc(tbl.contingency, plotit=FALSE))


###################################################
### code chunk number 12: explicateAUPC
###################################################
last.row <- nrow(tbl.contingency)
mid.row <- as.integer(round(last.row)/2)
print(tbl.contingency[c(1,mid.row,last.row),])


###################################################
### code chunk number 13: headtblInferred
###################################################
print(head(tbl.inferred, n=10))


###################################################
### code chunk number 14: g10targets
###################################################
print(subset(tbl.inferred, Regulator=="G10"))


###################################################
### code chunk number 15: g10inGoldStandard
###################################################
mtx.goldStandard["G10",]


###################################################
### code chunk number 16: simulateG10binding
###################################################
set.seed(37)
tbl.priors <- matrix(data=runif(100, 0, 0.4), nrow=10, ncol=10, 
                     dimnames=list(rownames(mtx.goldStandard), 
                                   colnames(mtx.goldStandard)))
tbl.priors["G10", c("G3", "G4")] <- runif(2, 0.8, 1.0)


###################################################
### code chunk number 17: inferWithPriors
###################################################
tbl.inferredWithPriors <- networkBMA(mtx.ts1, nTimePoints=nrow(mtx.ts1), 
                                     prior.prob=tbl.priors, self=FALSE)
tbl.inferredWithPriors <- 
   tbl.inferredWithPriors[order(tbl.inferredWithPriors$PostProb, 
                                decreasing=TRUE),]

print(tbl.inferredWithPriors)


###################################################
### code chunk number 18: aupr2
###################################################
tbl.contingencyWithPriors <- contabs.netwBMA(tbl.inferredWithPriors, tbl.goldStandard[,-3])
#precisionWithPriors <- with(tbl.contingencyWithPriors, TP/(TP+FP))
#recallWithPriors <- with(tbl.contingencyWithPriors, TP/(TP+FN))
plot(pr$recall, pr$precision, type='b',
     xlab='RECALL', ylab='PRECISION', 
     col=colors[1],
     xlim=c(0,1), ylim=c(0,1))

prWithPriors <- scores( tbl.contingencyWithPriors, what = c("precision","recall"))
lines(prWithPriors$recall, prWithPriors$precision, type='b',
      xlab='RECALL', ylab='PRECISION', 
      col=colors[2],
      xlim=c(0,1), ylim=c(0,1))

legend.titles = c("expression only", "with priors")
legend (0.6, 0.9, legend.titles, colors)

print(prc(tbl.contingencyWithPriors, plotit=FALSE))


###################################################
### code chunk number 19: printx
###################################################
print(tbl.contingencyWithPriors)


