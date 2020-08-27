##============================================================================##
##                                SeqGate Package                             ##
##============================================================================##

##  Copyright 2020 Christelle Reynès, Stéphanie Rialle
##  This file is part of SeqGate.
##  SeqGate is free software: you can redistribute it and/or modify it under
##  the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##  SeqGate is distributed in the hope that it will be useful, but WITHOUT ANY
##  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
##  FOR A PARTICULAR PURPOSE.
##  See the GNU General Public License for more details.
##  You should have received a copy of the GNU General Public License along with
##  this program.
##  If not, see <http://www.gnu.org/licenses/>.


## -----------------------------------------------------------------------------
## applySeqGate function
## Implements the SeqGate methods to filter lowly expressed genes
## -----------------------------------------------------------------------------

applySeqGate <- function(readCounts, conditions, prop0=NULL, percentile=NULL,
    propUpThresh=0.9) {
    
    if (missing(readCounts)) {
        stop("Missing mandatory argument readCounts.")
    }
    if (missing(conditions)) {
        stop("Missing mandatory argument conditions.")
    }    
    readCounts <- checkReadCounts(readCounts, conditions)
    prop0 <- checkAndDefineProp0(prop0, conditions)
    percentile <- checkAndDefinePercentile(percentile,conditions)
    checkPropUpThresh(propUpThresh)
    checkForReplicates(conditions)
    cat(getMessageWithParameters(prop0, percentile, propUpThresh))

    ## Getting the distribution of 'max' counts and computing threshold
    vecMax <- getDistMaxCounts(readCounts,conditions,prop0)
    if(length(vecMax)>0){
        threshold <- round(quantile(vecMax, percentile))
        cat(paste("The applied threshold is ", threshold, ".\n", sep=""))
    }else{
        stop("Impossible to calculate the distribution of max counts.")
    }

    ## Applying filter
    uniqConds <- unique(conditions)
    nbUpThresh <- rep(NA, length(uniqConds))
    for (i in seq_len(length(uniqConds))) { 
        nbUpThresh[i] <- ceiling(sum(conditions == uniqConds[i]) * propUpThresh)
    }
    applyFilter <- function(readCountsLoc, nbUpLoc = nbUpThresh,
        conditionsLoc = conditions, threshLoc = threshold) {
        countLoc <- tapply(readCountsLoc, conditionsLoc,
            function(vec,thresh)sum(vec >= thresh), thresh=threshLoc)
        if (any(countLoc >= nbUpLoc)) return(TRUE) else return(FALSE)
    }
    filterRes<-apply(readCounts, 1, applyFilter)
    
    ## Return the matrix of filtered and unfiltered genes 
    return(list(matrixKept = readCounts[filterRes,],
        matrixOut = readCounts[filterRes == FALSE,],
        onFilter = as.numeric(filterRes),
        threshold =  threshold))
}

## -----------------------------------------------------------------------------
## Unexported internal functions
## -----------------------------------------------------------------------------

## Checking readCounts type and content
checkReadCounts <- function(counts, conditions) {
    if(!is.matrix(counts)) {
        if(is.data.frame(counts)){
            isNum <- vapply(counts, is.numeric, FUN.VALUE=logical(1))
            i <- 1
            while(i <= length(isNum) & isNum[i]){
                i <- i+1
            }
            if(i > length(isNum)){
                counts <- as.matrix(counts)
                message("readCounts dataframe has been converted to a matrix.")
            }else{
                stop("readCounts contains non numeric columns. Please provide a 
matrix of read counts.")
            }
        }else{
            stop("readCounts is not a matrix.")
        }
    }else{
        if(!is.numeric(counts)){
            stop("readCounts contains non numeric data. Please provide a matrix 
of read counts.")
        }
    }
    
    if (ncol(counts) != length(conditions)) {
        stop("The number of columns in the readCounts matrix must be equals to 
the number of elements in the conditions vector.")
    }
    
    return(counts)
}

## Ckecking other arguments
checkAndDefineProp0 <- function(prop0, conditions) {
    if (is.null(prop0)) {
        prop0 <- (min(table(conditions))-1)/(min(table(conditions)))
    } else {
        if (prop0 <= 0 | prop0 >= 1) {
            stop("prop0 must be >0 and <1.")
        }
    }
    return(prop0)
}

checkAndDefinePercentile <- function(percentile, conditions) {
    if (is.null(percentile)) {
        if (min(table(conditions)) < 5) {
            percentile <- 0.9 
        } else {
            percentile <- max(c(0.5,-0.06*min(table(conditions))+1.1))
        }
    } else {
        if (percentile < 0 | percentile > 1) {
            stop("percentile must be comprised between 0 and 1.")
        } else {
            if (percentile == 0) {
                warning(
"percentile == 0, the filering threshold will probably be 0, meaning that no 
filtering will be realised."
                )
            }
        }
    }
    return(percentile)
}

checkPropUpThresh <- function(propUpThresh) {
    if (propUpThresh <= 0 | propUpThresh > 1) {
        stop("propUpThresh must be >0 and <=1.")
    }
}

## Checking for presence of replicates
checkForReplicates <- function(conditions){
    if(!any(duplicated(conditions))){
        stop("Condition vector does not contain replicates.")
    }
}

## Checking that the data does contain zeros
checkForZeros <- function(counts) {
    nTot <- nrow(counts) * ncol(counts)
    if (table(apply(counts,1,function(vec) vec == 0))["FALSE"] == nTot) {
        stop("Data does not contain zero counts. Sorry, SeqGate can't be 
applied.")
    }
}

## getMessageWithParameters
getMessageWithParameters <- function(prop0, percentile, propUpThresh) {
    return(paste("The SeqGate filter parameters applied are:\nprop0=", prop0, 
    "\npercentile=", percentile, "\npropUpThresh=",propUpThresh,"\n", sep = ""))
}

## Getting the distribution of 'max' counts
getDistMaxCounts <- function(readCounts, conditions, prop0){
    checkForZeros(readCounts)
    cpt0 <- apply(readCounts,1,sum)
    vecMax <- NULL
    for (i in unique(conditions)) {
        nb0Min <- ceiling(sum(conditions == i) * prop0)
        iC <- which(conditions == i)
        nb0 <- apply(readCounts[cpt0 != 0, iC], 1, function(vec) sum(vec == 0))
        iLow <- which(nb0 >= nb0Min)
        subReadCounts <- readCounts[which(cpt0 != 0)[iLow],iC]
        if(is.matrix(subReadCounts)){
            vecMax <- c(vecMax, 
                        apply(readCounts[which(cpt0 != 0)[iLow],iC], 1, max))
        }else{
            stop("Too few zeros along with counts in condition ",i,".")
        }
    }
    return(vecMax)
}
