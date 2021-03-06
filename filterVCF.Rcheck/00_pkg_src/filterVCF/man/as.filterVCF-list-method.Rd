% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_class_methods.R
\docType{methods}
\name{as.filterVCF,list-method}
\alias{as.filterVCF,list-method}
\title{An S4 methods for object list}
\usage{
\S4method{as.filterVCF}{list}(object)
}
\description{
Coerce a list to filterVCF object. Some slots can be empty.
}
\examples{
list_vcf <- list(
vcf = new("vcfR"),
vcf_filtered = new("vcfR"),
vcf_inv = GRanges(),
vcf_inv_filtered = GRanges()
)
new_filterVCF <- as.filterVCF(list_vcf)
}
