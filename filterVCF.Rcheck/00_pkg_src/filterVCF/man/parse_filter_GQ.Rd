% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_optionParser.R
\name{parse_filter_GQ}
\alias{parse_filter_GQ}
\title{parse options for filter by genotype quality}
\usage{
parse_filter_GQ(opt_filter_GQ)
}
\arguments{
\item{opt_filter_depth}{Character string, see help of filterVCF.R pipeline for more info}
}
\value{
list of threshold for filter on GQ and RGQ (invariant)
}
\description{
Intented to be use with the pileline filterVCF.R
}
\examples{

Parse the option from a call to 
option_list <- initialise_option()
opt <- parse_args(OptionParser(option_list=option_list))
opt$filter_GQ <- "5,10"
#OR
opt$filter_GQ <- "auto"
parse_filter_GQ(opt$filter_GQ)

}
\references{

}
\author{
~~Benjamin Laenen~~
}
\keyword{parser}
