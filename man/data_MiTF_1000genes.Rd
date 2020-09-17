\name{data_MiTF_1000genes}
\alias{data_MiTF_1000genes}

\docType{data}

\title{
    Extract of a Transciptome Dataset from Human Melanoma Cell Line
}

\description{
    This data set is an extract from a dataset from human transciptome studied 
    by Strub et al. (2011).
}

\usage{data(data_MiTF_1000genes)}
\format{
    A data frame with 1,000 raw gene counts on the following 6 samples.
    \describe{
        \item{\code{A1}}{a numeric vector}: counts from sample A1 (condition A)
        \item{\code{A2}}{a numeric vector}: counts from sample A2 (condition A)
        \item{\code{B1}}{a numeric vector}: counts from sample B1 (condition B)
        \item{\code{B2}}{a numeric vector}: counts from sample B2 (condition B)
        \item{\code{A3}}{a numeric vector}: counts from sample A3 (condition A)
        \item{\code{B3}}{a numeric vector}: counts from sample B3 (condition B)
    }
}
\details{
    The extract counts 1,000 genes measured in 6 samples, from 2 biological 
    conditions (cells expressing the Microphtalmia Transcription Factor (MiTF) 
    are compared to cells in which the MiTF is repressed). The raw read counts 
    have been retrieved from the Supplementary Materials of Dillies et al. 
    (2012). The initial dataset counts 37,720 genes.
}
\source{
    Dillies, M.A. et al. (2012) A comprehensive evaluation of normalization 
    methods for Illumina high-throughput RNA sequencing data analysis. Brief. 
    Bioinformatics 2013 Nov;14(6):671-83
}
\references{
    Strub, T., Giuliano, S., Ye T,, et al. Essential role of microphthalmia 
    transcription factor for DNA replication, mitosis and genomic stability in 
    melanoma. Oncogene 2011;30:2319--32
}
\examples{
    data(data_MiTF_1000genes)
}
