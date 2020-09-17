\name{applySeqGate}
\alias{applySeqGate}
\alias{seqGate}
\alias{seqgate}

\title{Performs a Filtering of Lowly Expressed Features}

\description{
    Implements the SeqGate methods to filter lowly expressed features (e.g. 
    genes)
}

\usage{
    applySeqGate(readCountsSE, assayName, colCond, prop0=NULL, percentile=NULL, 
                propUpThresh=0.9)
}

\arguments{
    \item{readCountsSE}{
        \code{SummarizedExperiment} object including the \code{"assayName"} 
        assay, that is a matrix of read counts. Rows corresponds to features 
        (genes or other genomic features) and columns corresponds to sample 
        libraries. 
        The read counts must contain zeros, as the SeqGate method is based on 
        the distribution of counts in replicates along with zero counts.
    }
    \item{assayName}{
        character string giving the name of the assay in \code{readCountsSE} for
        which we want to perform the filtering.
    }
    \item{colCond}{
        character string giving the name of the column describing the biological
        conditions to which the samples belong to, in the colData DataFrame 
        of \code{readCountsSE}.
    }
    \item{prop0}{
        minimal proportion of zeros among a condition to consider that the 
        feature is not or lowly expressed. \code{prop0} must be > 0 and < 1.
    }
    \item{percentile}{
        percentile used on the 'max' distribution to determine the filtering 
        threshold value. \code{percentile} must be >= 0 and <=1.
    }
    \item{propUpThresh}{
        proportion of counts to be above the threshold in at least one condition
        to keep the feature. \code{propUpThresh} must be > 0 and <= 1.
    }
}

\details{
    In order to find a theshold value to filter lowly expressed features, 
    SeqGate analyzes the distribution of counts found in replicates along with
    zero counts. More specifically, features with a proportion of at least 
    \code{prop0} zeros in one condition are selected. The distribution of counts
    found in replicates of that same condition along with those zeros is 
    computed. The chosen threshold is the count value corresponding to the 
    \code{percentile} of this distribution. Finally, features having a  
    proportion \code{propUpThresh} of replicates with counts below that value 
    in all conditions are filtered.
    
    If \code{prop0} is not provided, it is set to the number of replicates minus
    one divided by the max total number of replicates. For example, if one of
    the condition of the experiment counts 2 replicates and the other condition 
    4 replicates, the proportion \code{prop0} will be set to (4-1)/4=0.75.
    
    If \code{percentile} is not provided, it is set to 0.90 if at least one of 
    the conditions have 5 replicates or less. Otherwise, it is set to a value 
    comprised between 0.5 and 0.9, depending on the number or replicates.
}

\value{
    The input {\code{SummarizedExperiment}} with the following added elements{
    \item{rowData(readCountsSE)$onFilter}{
        vector which \code{length = nrow(readCounts)} indicating if features are
        kept after filtering (\code{value = TRUE}) or not (\code{value = FALSE})
    }
    \item{metadata(readCountsSE)$threshold}{
        applied filter threshold.
    }
}

\references{
Rialle, R. et al. (2020): SeqGate: a bioconductor package to perform data-driven
filtering of RNAseq datasets (manuscript in preparation)
}

\author{
    Christelle Reynès \email{christelle.reynes@igf.cnrs.fr},\cr
    Stéphanie Rialle \email{stephanie.rialle@mgx.cnrs.fr}
}

\examples{
    # Loading of input data frame
    data(data_MiTF_1000genes)
    # Annotating conditions
    cond<-c("A","A","B","B","A","B")
    # Setting the SummarizedExperiment input
    rowData <- DataFrame(row.names=rownames(data_MiTF_1000genes))
    colData <- DataFrame(Conditions=cond)
    counts_strub <- SummarizedExperiment(
                        assays=list(counts=data_MiTF_1000genes),
                        rowData=rowData,
                        colData=colData)
    # Applying SeqGate
    applySeqGate(counts_strub,"counts","Conditions")
    # Getting the matrix of kept genes after filtering
    keptGenes <- assay(counts_strub[counts_strub$onFilter == TRUE,])
}
