\name{SeqGate-package}

\alias{SeqGate-package}
\alias{SeqGate}

\docType{package}

\title{
\packageTitle{SeqGate}
}

\description{
    SeqGate is a method to filter lowly expressed features (e.g. genes).
}
\details{
    From a matrix of counts where lines correspond to features and columns to
    biological samples, provided as a SummarizedExperiment object, a threshold 
    is computed and applied in order to filter lowly expressed features. The 
    threshold is computed based on the distribution of counts measured along 
    with zeros within replicates of the same condition. The objective of SeqGate
    is to rationalize the filtering step by using the information of replicate 
    samples. The computed threshold corresponds to the count value below which 
    we can not be sure that the count can be considered different from zero.\cr
    The filtering is made by calling the \code{applySeqGate()} function.
}
\author{
Christelle Reynès \email{christelle.reynes@igf.cnrs.fr},\cr
Stéphanie Rialle \email{stephanie.rialle@mgx.cnrs.fr},\cr
Maintainer: \packageMaintainer{SeqGate}
}
\references{
Rialle, R. et al. (2020): SeqGate: a bioconductor package to perform data-driven
filtering of RNAseq datasets \emph{manuscript in preparation}
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
