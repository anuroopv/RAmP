% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.filter.R
\name{QC.filter}
\alias{QC.filter}
\title{QC.filter}
\usage{
QC.filter(data, filter.thr = NA, filter.protein.min = NULL)
}
\arguments{
\item{data}{output file (or output) from editData function}

\item{filter.thr}{Only if filter.protein.type = condition. Numerical value less than the number of relicates in either condition (Ex. 0 indicates the protein should have no NAs in all replicates of atleast one condition while 1 indicates they can have one NAs)}

\item{filter.protein.min}{Only if filter.protein.type = fraction. Any value between 0-1 Any value between 0-1 (Ex. 0.75 indicates the protein should not have NAs in 75\% of all samples)}
}
\value{
Filtered and vsn normalized data
}
\description{
Performs quality control and filters the missing values based on parameters provided
}
\examples{
norm_data <- QC.filter(data = lfqdata, filter.protein.type = "fraction", filter.protein.min = 0.75)
}
