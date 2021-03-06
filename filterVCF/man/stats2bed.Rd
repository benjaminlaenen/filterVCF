% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_popstat.R
\name{stats2bed}
\alias{stats2bed}
\title{Save output from Popstats into a BED file}
\usage{
stats2bed(gr, outputdir = ".", output_name = "stat.bed")
}
\arguments{
\item{gr}{a GRange or output from Popstat()}

\item{outputdir="."}{Output directory}

\item{output_name="stat.bed"}{name of the file}
}
\value{
Save a BEDfile
}
\description{
More generally this function saves a Granges with metadata to a GRange and
is equivalent to Grange2bed(gr, output_name = "stat.bed", merge = FALSE,
keep_extra_col = TRUE)
}
\examples{


stats2bed(Popstat(filterVCFobject), windows_size = 20000)

}
\references{

}
\author{
~~Benjamin Laenen~~
}
