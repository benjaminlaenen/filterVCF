% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_class_methods.R
\docType{package}
\name{filterVCF-package}
\alias{filterVCF}
\alias{filterVCF-package}
\title{filterVCF: filterVCF}
\description{
Filter Variant Call Format (VCF) file and compute population genomic
}
\details{
This package is intented to filter large VCF file using a panel of filter
 based on Metadata, coverage, genotype, genotype quality, repeats. It introduce
 two new S4 classes "filterVCF" and "VIF" (VariantInvariantFixed) that keep
 track of the all the sites that were sequenced in a condensed and memmory
 efficient way using Genomic Ranges and ffdf dataframe for large object. The
 rational behind the package is that not all sites in the genome are available
 for analysis, and it is difficult to keep track of which sites were filtered
 or not and how many remains after filtering. The package always ties SNPs,
 invariant and fixed sites together and offer many utilities to manipulates the
 two type of object. Intersection, substraction, splitting, subsampling and
 computation of population genetics statistics are supported. The filtering
 applied consist of hard filter, such as depth cutoff. Modern tools such as
 GATK HaplotypeCaller use machine learning algorithm to create custom filter
 based on an curated set of known SNPs. This methods should be more robust and
 should become standart in the future, however for non-model or model species
 with less extensive knowledge on true variant or complex genomes (e.g. higly
 repetitive plant genome) hard filtering is still a prerequisite before doing
 any population genetic analysis.The package integrates with many other
 libraries of BioConductor and do not try to re-invent the wheel but use widely
 used object of general class such as GRanges, DNAStringSet and vcfR. It thus
 integrates easilly in a BioConductor workflow. The package is also used in a
 pipeline called filterVCF.R that automates the filtering of a raw VCF file
 containing all sites (tested on GATK VCF file) to be run on high performance
 cluster in UNIX. This pipeline handles all the steps of reading, filtering,
 saving output and create HTML reports.
}
\keyword{internal}
