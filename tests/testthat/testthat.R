context("check applySeqGate function")
library(SeqGate)

## Load test data

load("dataDillies100.Rda")
cond<-c("A","A","B","B","A","B")
load("dataNo0.Rda")
load("dataDec.Rda")
badcond<-c("A","B","C","D","E","F")

readCounts<-matrix(c(rep(c(0,10,20),4)),nrow=3,ncol=4)

## Perfom tests

test_that("applySeqGate function gives appropriate results in normal usage", {
    expect_output(res <- applySeqGate(dataDillies,cond),"The SeqGate filter parameters applied are")
    expect_output(resD <- applySeqGate(dataDec,cond),"The SeqGate filter parameters applied are")
    expect_is(res,"list")
    expect_output(applySeqGate(dataDillies,cond,0.5,0.4,0.6),"The SeqGate filter parameters applied are")
    expect_is(res$matrixKept,"matrix")
    expect_is(res$matrixOut,"matrix")
    expect_is(res$onFilter,"numeric")
    expect_equal(nrow(res$matrixKept)+nrow(res$matrixOut),nrow(dataDillies))
    expect_equal(intersect(rownames(res$matrixOut),rownames(res$matrixKept)),character(0))
    expect_equal(sum(res$onFilter),nrow(res$matrixKept))
    expect_equal(length(res$onFilter),(nrow(res$matrixKept)+nrow(res$matrixOut)))
})

test_that("applySeqGate function gives appropriate results in extreme usage", {
    expect_error(applySeqGate(dataNo0,cond),"Data does not contain zero counts")
    expect_error(applySeqGate(matrix(c(rep(c(0,10,20),6)),nrow=3,ncol=6),cond),"Impossible to calculate the distribution of max counts.")
})

test_that("applySeqGate function returns warning in extreme or weird situations", {
    expect_output(applySeqGate(dataDillies,cond,percentile=0),"The SeqGate filter parameters applied are")
    expect_warning(applySeqGate(dataDillies,cond,percentile=0),"percentile == 0")
    expect_output(applySeqGate(dataDillies,cond,percentile=1),"The SeqGate filter parameters applied are")
    expect_output(applySeqGate(dataDillies,cond,propUpThresh=1),"The SeqGate filter parameters applied are")
})

test_that("applySeqGate returns the good error messages", {
    data2 <- cbind.data.frame(rownames(dataDillies),dataDillies)
    expect_error(applySeqGate(data2,cond),"readCounts contains non numeric columns")
    expect_error(applySeqGate(dataDillies,c("A","A","B","B","A")),"number of elements in the conditions vector")
    expect_error(applySeqGate(dataDillies),"Missing mandatory argument")
    expect_error(applySeqGate(cond),"Missing mandatory argument")
    expect_error(applySeqGate(),"Missing mandatory argument")
    expect_error(applySeqGate(dataDillies,cond,prop0=50),"prop0 must be >0 and <1")
    expect_error(applySeqGate(dataDillies,cond,prop0=0),"prop0 must be >0 and <1")
    expect_error(applySeqGate(dataDillies,cond,prop0=1),"prop0 must be >0 and <1")
    expect_error(applySeqGate(dataDillies,cond,percentile=95),"percentile must be comprised between 0 and 1.")
    expect_error(applySeqGate(dataDillies,cond,propUpThresh=0),"propUpThresh must be >0 and <=1.")
    expect_error(applySeqGate(dataDillies,cond,propUpThresh=-5),"propUpThresh must be >0 and <=1.")
    expect_error(applySeqGate(dataDillies,badcond),"Condition vector does not contain replicates.")
})
