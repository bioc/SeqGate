context("check applySeqGate function")
library(SeqGate)

## Generating test data

cond<-c("A","A","B","B","A","B")
nrow=100
ncol=6
counts <- rnbinom(n = (nrow*ncol), size = 2, mu=100)
counts_with_zeros<-matrix(counts,nrow=nrow)
for(i in seq(1,nrow,round(nrow/10))){
    counts_with_zeros[i,]<-c(0,0,100,150,i,80)
    counts_with_zeros[i,]<-c(50,60,(10+i),0,55,0)
}
counts_with_zeros[c(2,12),] <- rep(0,ncol)
rownames(counts_with_zeros) <- paste("gene",c(1:nrow),sep="")

counts_no_zeros <- matrix(c(rep(c(5,10,20),4)),nrow=3,ncol=ncol)

counts_dec <- matrix(c(rep(c(5.1,0,10.2,20.1,6.4,7.9,40.5),6)),nrow=7,ncol=ncol)
counts_dec[1,c(2,5)]<-c(0,0)
counts_dec[5,c(1,2)]<-c(0,0)
counts_dec[3,c(4,6)]<-c(0,0)
counts_dec[6,c(3,6)]<-c(0,0)

counts_few_zeros <- matrix(c(rep(c(5,0,10,20),6)),nrow=4,ncol=ncol)
counts_few_zeros[1,c(2,5)]<-c(0,0)
counts_few_zeros[3,c(4,6)]<-c(0,0)

badcond<-c("A","B","C","D","E","F")

## Perfom tests

test_that("applySeqGate function gives appropriate results in normal usage", {
    expect_output(res <- applySeqGate(counts_with_zeros,cond),"The SeqGate filter parameters applied are")
    expect_output(resD <- applySeqGate(counts_dec,cond),"The SeqGate filter parameters applied are")
    expect_is(res,"list")
    expect_output(applySeqGate(counts_with_zeros,cond,0.5,0.4,0.6),"The SeqGate filter parameters applied are")
    expect_is(res$matrixKept,"matrix")
    expect_is(res$matrixOut,"matrix")
    expect_is(res$onFilter,"numeric")
    expect_equal(nrow(res$matrixKept)+nrow(res$matrixOut),nrow(counts_with_zeros))
    expect_equal(intersect(rownames(res$matrixOut),rownames(res$matrixKept)),character(0))
    expect_equal(sum(res$onFilter),nrow(res$matrixKept))
    expect_equal(length(res$onFilter),(nrow(res$matrixKept)+nrow(res$matrixOut)))
})

test_that("applySeqGate function gives appropriate results in extreme usage", {
    expect_error(applySeqGate(counts_no_zeros,cond),"Data does not contain zero counts")
    expect_error(applySeqGate(matrix(c(rep(c(0,10,20),6)),nrow=3,ncol=6),cond),"Impossible to calculate the distribution of max counts.")
})

test_that("applySeqGate function returns warning in extreme or weird situations", {
    expect_warning(applySeqGate(counts_with_zeros,cond,percentile=0),"percentile == 0")
    expect_output(applySeqGate(counts_with_zeros,cond,percentile=1),"The SeqGate filter parameters applied are")
    expect_output(applySeqGate(counts_with_zeros,cond,propUpThresh=1),"The SeqGate filter parameters applied are")
})

test_that("applySeqGate returns the good error messages", {
    data2 <- cbind.data.frame(rownames(counts_with_zeros),counts_with_zeros)
    expect_error(applySeqGate(data2,cond),"readCounts contains non numeric columns")
    expect_error(applySeqGate(counts_with_zeros,c("A","A","B","B","A")),"number of elements in the conditions vector")
    expect_error(applySeqGate(counts_with_zeros),"Missing mandatory argument")
    expect_error(applySeqGate(cond),"Missing mandatory argument")
    expect_error(applySeqGate(),"Missing mandatory argument")
    expect_error(applySeqGate(counts_with_zeros,cond,prop0=50),"prop0 must be >0 and <1")
    expect_error(applySeqGate(counts_with_zeros,cond,prop0=0),"prop0 must be >0 and <1")
    expect_error(applySeqGate(counts_with_zeros,cond,prop0=1),"prop0 must be >0 and <1")
    expect_error(applySeqGate(counts_with_zeros,cond,percentile=95),"percentile must be comprised between 0 and 1.")
    expect_error(applySeqGate(counts_with_zeros,cond,propUpThresh=0),"propUpThresh must be >0 and <=1.")
    expect_error(applySeqGate(counts_with_zeros,cond,propUpThresh=-5),"propUpThresh must be >0 and <=1.")
    expect_error(applySeqGate(counts_with_zeros,badcond),"Condition vector does not contain replicates.")
    expect_error(applySeqGate(counts_few_zeros,cond),"Too few zeros along with counts in condition")
})
