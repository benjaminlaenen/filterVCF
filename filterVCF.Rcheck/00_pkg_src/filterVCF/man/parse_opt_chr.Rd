% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_optionParser.R
\name{parse_opt_chr}
\alias{parse_opt_chr}
\title{parse options for selecting chromosome}
\usage{
parse_opt_chr(opt_chr)
}
\arguments{
\item{opt_chr}{Character string, see help of filterVCF.R pipeline for more info}
}
\value{
Character string with chromosome name
}
\description{
Intented to be use with the pileline filterVCF.R
}
\examples{

Parse the option from a call to 
option_list <- initialise_option()
opt <- parse_args(OptionParser(option_list=option_list))
opt$chr <- "chr1,chr2, chr3"
parse_opt_chr(opt$chr)

}
\references{

}
\author{
~~Benjamin Laenen~~
}
\keyword{parser}
