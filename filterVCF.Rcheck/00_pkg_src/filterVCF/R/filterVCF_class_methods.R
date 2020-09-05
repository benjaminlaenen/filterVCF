# ==========================================================================
#     Package description
# ==========================================================================
#' @details
#'  This package is intented to filter large VCF file using a panel of filter
#'  based on Metadata, coverage, genotype, genotype quality, repeats. It introduce
#'  two new S4 classes "filterVCF" and "VIF" (VariantInvariantFixed) that keep
#'  track of the all the sites that were sequenced in a condensed and memmory
#'  efficient way using Genomic Ranges and ffdf dataframe for large object. The
#'  rational behind the package is that not all sites in the genome are available
#'  for analysis, and it is difficult to keep track of which sites were filtered
#'  or not and how many remains after filtering. The package always ties SNPs,
#'  invariant and fixed sites together and offer many utilities to manipulates the
#'  two type of object. Intersection, substraction, splitting, subsampling and
#'  computation of population genetics statistics are supported. The filtering
#'  applied consist of hard filter, such as depth cutoff. Modern tools such as
#'  GATK HaplotypeCaller use machine learning algorithm to create custom filter
#'  based on an curated set of known SNPs. This methods should be more robust and
#'  should become standart in the future, however for non-model or model species
#'  with less extensive knowledge on true variant or complex genomes (e.g. higly
#'  repetitive plant genome) hard filtering is still a prerequisite before doing
#'  any population genetic analysis.The package integrates with many other
#'  libraries of BioConductor and do not try to re-invent the wheel but use widely
#'  used object of general class such as GRanges, DNAStringSet and vcfR. It thus
#'  integrates easilly in a BioConductor workflow. The package is also used in a
#'  pipeline called filterVCF.R that automates the filtering of a raw VCF file
#'  containing all sites (tested on GATK VCF file) to be run on high performance
#'  cluster in UNIX. This pipeline handles all the steps of reading, filtering,
#'  saving output and create HTML reports.
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"
# ==========================================================================
#                    FUNCTIONS
# ==========================================================================

# ==========================================================================
#                SET Class filterVCFobject
# ==========================================================================
# Set a virtual class to be able
# to inherit methods between filterVCF and VIF
#' @export
setClass(Class="ParentClass", representation=representation("VIRTUAL"))

#' An S4 class to filterVCF object.
#'
#' @slot vcf The original vcf with only snp of class "vcfR",
#' @slot vcf_filtered The filtered vcf of class "vcfR",
#' @slot vcf_inv The list of invariant sites of class "GRanges",
#' @slot vcf_inv_filtered The list of invariant filtered sites of class "GRanges",
#' @slot fixed_inv_Grange The list of invariant fixed sites of class"GRanges",
#' @slot reference A DNAStringSet object with the reference used to create the object,
#' @slot master_filter All sites too filter for snps. A logical vector with the same length as the @vcf,
#' @slot master_filter_inv  All sites too filter for invariant. A logical vector with the same length as the @vcf_inv,
#' @slot filters A list of all the filters applied to snp. List of logical vector with length equal to @vcf,
#' @slot filters_inv A list of all the filters applied to invariants. List of logical vector with length equal to @vcf,
#' @slot removed_sites All sites removed in snp of class "GRanges",
#' @slot removed_sites_inv All sites removed in invariant of class "GRanges",
#' @slot bed_GRange_merged Merged sites from the bedfiles inputed for filtering of class "GRanges,
#' @slot windows_with_half_repeats_to_remove Windows with 50 percent of repeats of class "GRanges",
#' @slot fix_het_pop Sites that are fixed heterozygotes in the selected population of class "GRanges",
#' @slot filter_DP_repeats  A list of filter based on local standardised depth. A list of logical vector of length @vcf,
#' @slot stats_vcf_filtered Unused for the moment. Class "character",
#' @slot opt List of option used by the pipeline filterVCF.R,
#' @export
#' @import IRanges GenomicRanges rtracklayer Biostrings
setClass("filterVCF", slots = c(
vcf                                 = "vcfR",
vcf_filtered                        = "vcfR",
vcf_inv                             = "GRanges",
vcf_inv_filtered                    = "GRanges",
fixed_inv_Grange                    = "GRanges",
reference                           = "DNAStringSet",
master_filter                       = "logical",
master_filter_inv                   = "logical",
filters                             = "list",
filters_inv                         = "list",
removed_sites                       = "GRanges",
removed_sites_inv                   = "GRanges",
bed_GRange_merged                   = "GRanges",
windows_with_half_repeats_to_remove = "GRanges",
fix_het_pop                         = "GRanges",
filter_DP_repeats                   = "list",
stats_vcf_filtered                  = "character",
opt 								= "list"),
contains = "ParentClass"
)


# ==========================================================================
#                            Generic
# ==========================================================================

setGeneric("show.filterVCF", function(object) standardGeneric("show.filterVCF"))
setGeneric("as.filterVCF", function(object) standardGeneric("as.filterVCF"))
setGeneric("as.list.filterVCF", function(object) standardGeneric("as.list.filterVCF"))
setGeneric("repair", function(object) standardGeneric("repair"))


# ==========================================================================
#                           Methods
# ==========================================================================

#' An S4 methods for object filterVCF
#'
#' Print information on the current object.
#' @export
setMethod("show", signature(object = "filterVCF"), function(object) {
	nsample <- ncol(vcfRaw(object)@gt) - 1
	nsample_filtered <- ncol(vcf(object)@gt) - 1
	if(nsample != nsample_filtered) message(paste(strwrap(sprintf("Warnings : Number of sample is not the same for filtered VCF (%s) and original VCF (%s). This is normal if output were merged after running chunks. Original VCF should thus not be use but just kept for the position.\n", nsample_filtered, nsample), width = 120,initial = "\n\t", prefix="\n\t"), collapse =""))
	
	n_SNP <- nrow(vcfRaw(object))
	n_SNP_filtered <- nrow(vcf(object))
	n_removed <- sum(masterFilter(object))	
	#stopifnot(n_SNP - n_SNP_filtered == n_removed)
	SNP <- c(raw = n_SNP, filtered = n_SNP_filtered, removed_sites = n_removed)

	n_inv <- sum(width(invariantRaw(object)))
	n_inv_filtered <- sum(width(invariant(object)))
	n_removed_inv <- sum(masterFilterInv(object))
	#stopifnot(n_inv - n_inv_filtered == n_removed_inv)
	inv <- c(raw = n_inv, filtered = n_inv_filtered, removed_sites = n_removed_inv)
	
	n_fixed <- sum(width(fixed(object)))
	names_filters <- names(Filters(object))
	n_per_filter <- as.integer(sapply(Filters(object), sum, USE.NAMES = FALSE))
	filterTable <- as.data.frame(cbind(Filter=names_filters, `Number of sites`=n_per_filter))

	printTable <- rbind(SNPs=SNP, Invariable =inv, Fixed = c(n_fixed, NA, NA))
	cat("
	# ==========================================================================
	#                         Object of class filterVCF
	# ==========================================================================\n")
	cat(paste("The number of sample is : ", nsample_filtered))
	cat("\n\n")
	cat("Number of sites before and after filtering")
	cat("\n\n")
	print(printTable)
	cat("List of the filters applied")
	cat("\n\n")
	print(filterTable)
	cat("\n")
	cat("
	# ==========================================================================
	#                        
	# ==========================================================================\n")
})


#' An S4 methods for object list
#'
#' Coerce a list to filterVCF object. Some slots can be empty.
#' @examples
#' list_vcf <- list(
#' vcf = new("vcfR"),
#' vcf_filtered = new("vcfR"),
#' vcf_inv = GRanges(),
#' vcf_inv_filtered = GRanges()
#' )
#' new_filterVCF <- as.filterVCF(list_vcf)
#' @export
setMethod("as.filterVCF", signature(object = "list"), function(object) {
	res <- filterVCF()
	res@vcf                                        = if(is.null(object$vcf)) new("vcfR") else object$vcf
	res@vcf_filtered                               = if(is.null(object$vcf_filtered)) new("vcfR") else object$vcf_filtered
	res@vcf_inv                                    = if(is.null(object$vcf_inv)) GRanges() else object$vcf_inv
	res@vcf_inv_filtered                           = if(is.null(object$vcf_inv_filtered)) GRanges() else object$vcf_inv_filtered
	res@fixed_inv_Grange                           = if(is.null(object$fixed_inv_Grange)) GRanges() else object$fixed_inv_Grange
	res@stats_vcf_filtered                         = as.character(object$stats_vcf_filtered)
	res@master_filter                              = if(is.null(object$master_filter)) logical()  else object$master_filter
	res@master_filter_inv                          = if(is.null(object$master_filter_inv)) logical()  else object$master_filter_inv
	res@filters                                    = if(is.null(object$filters)) list()  else object$filters
	res@filters_inv                                = if(is.null(object$filters_inv))  list() else object$filters_inv
	res@removed_sites                              = if(is.null(object$removed_sites))  GRanges() else object$removed_sites
	res@removed_sites_inv                          = if(is.null(object$removed_sites_inv)) GRanges()  else object$removed_sites_inv
	res@bed_GRange_merged                          = if(is.null(object$bed_GRange_merged)) GRanges() else object$bed_GRange_merged
	res@windows_with_half_repeats_to_remove        = if(is.null(object$windows_with_half_repeats_to_remove))GRanges() else object$windows_with_half_repeats_to_remove
	res@fix_het_pop                                = if(is.null(object$fix_het_pop)) GRanges()else object$fix_het_pop
	res@filter_DP_repeats                          = if(is.null(object$filter_DP_repeats)) list(GRanges()) else object$filter_DP_repeats
	res@reference                                  = if(is.null(object$reference)) new("DNAStringSet") else object$reference
	res@opt                                        = if(is.null(object$opt)) parse_args(OptionParser(option_list=initialise_option())) else object$opt
	return(res)
})


#' An S4 methods for object filterVCF
#'
#' Coerce a filterVCF to list object. If slots are empty they are filled up by object of the right class.
#' @export
setMethod("as.list.filterVCF", signature(object = "filterVCF"), function(object) {
	res <- list()
	res$vcf                                        = if(.hasSlot(object, "vcf")) object@vcf else new("vcfR")
	res$vcf_filtered                               = if(.hasSlot(object, "vcf_filtered")) object@vcf_filtered else new("vcfR") 
	res$vcf_inv                                    = if(.hasSlot(object, "vcf_inv")) object@vcf_inv else new("GRanges")
	res$vcf_inv_filtered                           = if(.hasSlot(object, "vcf_inv_filtered")) object@vcf_inv_filtered else new("GRanges")
	res$fixed_inv_Grange                           = if(.hasSlot(object, "fixed_inv_Grange")) object@fixed_inv_Grange else new("GRanges")
	res$reference                                  = if(.hasSlot(object, "reference")) object@reference else new("DNAStringSet")
	res$master_filter                              = if(.hasSlot(object, "master_filter")) object@master_filter else new("logical")
	res$master_filter_inv                          = if(.hasSlot(object, "master_filter_inv")) object@master_filter_inv else new("logical")
	res$filters                                    = if(.hasSlot(object, "filters")) object@filters else new("list")
	res$filters_inv                                = if(.hasSlot(object, "filters_inv")) object@filters_inv else new("list")
	res$removed_sites                              = if(.hasSlot(object, "removed_sites")) object@removed_sites else new("GRanges")
	res$removed_sites_inv                          = if(.hasSlot(object, "removed_sites_inv")) object@removed_sites_inv else new("GRanges")
	res$bed_GRange_merged                          = if(.hasSlot(object, "bed_GRange_merged")) object@bed_GRange_merged else new("GRanges")
	res$windows_with_half_repeats_to_remove        = if(.hasSlot(object, "windows_with_half_repeats_to_remove")) object@windows_with_half_repeats_to_remove  else new("GRanges")
	res$fix_het_pop                                = if(.hasSlot(object, "fix_het_pop")) object@fix_het_pop else new("GRanges")
	res$filter_DP_repeats                          = if(.hasSlot(object, "filter_DP_repeats")) object@filter_DP_repeats else new("GRanges")
	res$opt                                        = if(.hasSlot(object, "opt")) object@opt else parse_args(OptionParser(option_list=initialise_option()))
	res$stats_vcf_filtered                         = if(.hasSlot(object, "stats_vcf_filtered")) object@stats_vcf_filtered else NA_character_
	return(res)
})


#' An S4 methods for object filterVCF
#'
#' Repair an object by filling the missing slot with empty object. Function for compatibility if the definition changes between version.
#' @export
setMethod("repair", signature(object = "filterVCF"), function(object) {
	return(as.filterVCF(as.list(object)))
})


#' Initialize a filterVCF
#'
#' Initialize a filterVCF
#' @export
setMethod("initialize", "filterVCF", 
    function(.Object, ...) {
      .Object <- callNextMethod()
      .Object@opt <- parse_args(OptionParser(option_list=initialise_option()), args = NULL)
      .Object
    }
)

# ==========================================================================
#                          Accessor
# ==========================================================================

# Acces the filtered vcf
setGeneric("vcf", function(x) standardGeneric("vcf"))
#' @export
setMethod("vcf", signature = "ParentClass", function(x) x@vcf_filtered)

setGeneric("vcf<-", function(x, value) standardGeneric("vcf<-"))
#' @export
setMethod("vcf<-", signature = "ParentClass", function(x, value) {
  x@vcf_filtered <- value
  validObject(x)
  x
})


# Acces the original vcf
setGeneric("vcfRaw", function(x) standardGeneric("vcfRaw"))
#' @export
setMethod("vcfRaw", signature = "filterVCF", function(x) x@vcf)

setGeneric("vcfRaw<-", function(x, value) standardGeneric("vcfRaw<-"))
#' @export
setMethod("vcfRaw<-", signature = "filterVCF", function(x, value) {
  x@vcf <- value
  validObject(x)
  x
})


# Acces the invariant
setGeneric("invariant", function(x) standardGeneric("invariant"))
#' @export
setMethod("invariant", signature = "ParentClass", function(x) x@vcf_inv_filtered)

setGeneric("invariant<-", function(x, value) standardGeneric("invariant<-"))
#' @export
setMethod("invariant<-", signature = "ParentClass", function(x, value) {
  x@vcf_inv_filtered <- value
  validObject(x)
  x
})

# Acces the original invariant
setGeneric("invariantRaw", function(x) standardGeneric("invariantRaw"))
#' @export
setMethod("invariantRaw", signature = "filterVCF", function(x) x@vcf_inv)

setGeneric("invariantRaw<-", function(x, value) standardGeneric("invariantRaw<-"))
#' @export
setMethod("invariantRaw<-", signature = "filterVCF", function(x, value) {
  x@vcf_inv <- value
  validObject(x)
  x
})

# Acces the fixed sites
setGeneric("fixed", function(x) standardGeneric("fixed"))
#' @export
setMethod("fixed", signature = "ParentClass", function(x) x@fixed_inv_Grange)

setGeneric("fixed<-", function(x, value) standardGeneric("fixed<-"))
#' @export
setMethod("fixed<-", signature = "ParentClass", function(x, value) {
  x@fixed_inv_Grange <- value
  validObject(x)
  x
})

# Acces the reference
setGeneric("reference", function(x) standardGeneric("reference"))
#' @export
setMethod("reference", signature = "ParentClass", function(x) x@reference)

setGeneric("reference<-", function(x, value) standardGeneric("reference<-"))
#' @export
setMethod("reference<-", signature = "ParentClass", function(x, value) {
  x@reference <- value
  validObject(x)
  x
})


# Access the master filter for snp (in regard to the original vcf). If TRUE the position do NOT pass the filter.
setGeneric("masterFilter", function(x) standardGeneric("masterFilter"))
#' @export
setMethod("masterFilter", signature = "filterVCF", function(x) x@master_filter)

setGeneric("masterFilter<-", function(x, value) standardGeneric("masterFilter<-"))
#' @export
setMethod("masterFilter<-", signature = "filterVCF", function(x, value) {
  x@master_filter <- value
  validObject(x)
  x
})

# Acces the master filter for invariant (in regard to the original vcf). If TRUE the position do NOT pass the filter.
setGeneric("masterFilterInv", function(x) standardGeneric("masterFilterInv"))
#' @export
setMethod("masterFilterInv", signature = "filterVCF", function(x) x@master_filter_inv)

setGeneric("masterFilterInv<-", function(x, value) standardGeneric("masterFilterInv<-"))
#' @export
setMethod("masterFilterInv<-", signature = "filterVCF", function(x, value) {
  x@master_filter_inv <- value
  validObject(x)
  x
})


# Acces the list of all Filters (in regard to the original vcf). If TRUE the position do NOT pass the filter.
setGeneric("Filters", function(x) standardGeneric("Filters"))
#' @export
setMethod("Filters", signature = "filterVCF", function(x) x@filters)

setGeneric("Filters<-", function(x, value) standardGeneric("Filters<-"))
#' @export
setMethod("Filters<-", signature = "filterVCF", function(x, value) {
  x@filters <- value
  validObject(x)
  x
})

# Acces the list of all Filters for invariant (in regard to the original vcf). If TRUE the position do NOT pass the filter.
setGeneric("FiltersInv", function(x) standardGeneric("FiltersInv"))
#' @export
setMethod("FiltersInv", signature = "filterVCF", function(x) x@filters_inv)

setGeneric("FiltersInv<-", function(x, value) standardGeneric("FiltersInv<-"))
#' @export
setMethod("FiltersInv<-", signature = "filterVCF", function(x, value) {
  x@filters_inv <- value
  validObject(x)
  x
})

# Acces the mask create by Bedfile(s) merged together
setGeneric("BEDFilter", function(x) standardGeneric("BEDFilter"))
#' @export
setMethod("BEDFilter", signature = "filterVCF", function(x) x@bed_GRange_merged)

setGeneric("BEDFilter<-", function(x, value) standardGeneric("BEDFilter<-"))
#' @export
setMethod("BEDFilter<-", signature = "filterVCF", function(x, value) {
  x@bed_GRange_merged <- value
  validObject(x)
  x
})



# Acces the removed sites for snp as GRange.
setGeneric("removed_sites", function(x) standardGeneric("removed_sites"))
#' @export
setMethod("removed_sites", signature = "ParentClass", function(x) x@removed_sites)

setGeneric("removed_sites<-", function(x, value) standardGeneric("removed_sites<-"))
#' @export
setMethod("removed_sites<-", signature = "ParentClass", function(x, value) {
  x@removed_sites <- value
  validObject(x)
  x
})

# Acces the removed sites for invariant as GRange.
setGeneric("removed_sites_inv", function(x) standardGeneric("removed_sites_inv"))
#' @export
setMethod("removed_sites_inv", signature = "ParentClass", function(x) x@removed_sites_inv)

setGeneric("removed_sites_inv<-", function(x, value) standardGeneric("removed_sites_inv<-"))
#' @export
setMethod("removed_sites_inv<-", signature = "ParentClass", function(x, value) {
  x@removed_sites_inv <- value
  validObject(x)
  x
})

# Acces the Grange object containing windows that do not pass the filter. 
# The filter is that half (hardcoded) of the windows interesect with the given Bedfile(e.g. repeats) 
setGeneric("Filter_WindowsProp", function(x) standardGeneric("Filter_WindowsProp"))
#' @export
setMethod("Filter_WindowsProp", signature = "ParentClass", function(x) x@windows_with_half_repeats_to_remove)

setGeneric("Filter_WindowsProp<-", function(x, value) standardGeneric("Filter_WindowsProp<-"))
#' @export
setMethod("Filter_WindowsProp<-", signature = "ParentClass", function(x, value) {
  x@windows_with_half_repeats_to_remove <- value
  validObject(x)
  x
})


# Acces the different filters for Normalized depth filtering accross sample and per sample rule in GRange format. The intervals are the one to remove.
setGeneric("Normalized_DP_Filter", function(x) standardGeneric("Normalized_DP_Filter"))
#' @export
setMethod("Normalized_DP_Filter", signature = "ParentClass", function(x) x@filter_DP_repeats)

setGeneric("Normalized_DP_Filter<-", function(x, value) standardGeneric("Normalized_DP_Filter<-"))
#' @export
setMethod("Normalized_DP_Filter<-", signature = "ParentClass", function(x, value) {
  x@filter_DP_repeats <- value
  validObject(x)
  x
})

# Acces windows with fixed heterozygous in a reference populations as a GRange object. The intervals are the one to remove.
setGeneric("FixHet_RefPop", function(x) standardGeneric("FixHet_RefPop"))
#' @export
setMethod("FixHet_RefPop", signature = "ParentClass", function(x) x@fix_het_pop)

setGeneric("FixHet_RefPop<-", function(x, value) standardGeneric("FixHet_RefPop<-"))
#' @export
setMethod("FixHet_RefPop<-", signature = "ParentClass", function(x, value) {
  x@fix_het_pop <- value
  validObject(x)
  x
})

# Access the options
setGeneric("GetOpt", function(x) standardGeneric("GetOpt"))
#' @export
setMethod("GetOpt", signature = "ParentClass", function(x) x@opt)

setGeneric("GetOpt<-", function(x, value) standardGeneric("GetOpt<-"))
#' @export
setMethod("GetOpt<-", signature = "ParentClass", function(x, value) {
  x@opt <- value
  validObject(x)
  x
})

# select one chromosome
setGeneric("CHROMselect", function(x, value) standardGeneric("CHROMselect"))
#' @export
setMethod("CHROMselect", signature = "ParentClass", function(x, value) {
  index_vcf <- getCHROM(vcf(x)) %in% value
  vcf(x) <- vcf(x)[index_vcf]

  invariant(x) <-  invariant(x)[seqnames(invariant(x)) == value]
  fixed(x) <-  fixed(x)[seqnames(fixed(x)) == value]
  removed_sites(x) <-  removed_sites(x)[seqnames(removed_sites(x)) == value]
  removed_sites_inv(x) <-  removed_sites_inv(x)[seqnames(removed_sites_inv(x)) == value]

  # For VIF object
  # specific slots
  if(.hasSlot(x, "divergence")){
	  divergence(x) <-  divergence(x)[seqnames(divergence(x)) == value]
	  repolKey(x) <-  repolKey(x)[seqnames(repolKey(x)) == value]
	  alternate(x) <-  alternate(x)[seqnames(alternate(x)) == value]
	  #reference
	  #boostrap
  } 
  if(.hasSlot(x, "master_filter")){
	  masterFilter(x) <- masterFilter(x)[index_vcf]
	  masterFilterInv(x) <- masterFilterInv(x)[as.logical(seqnames(invariant(x)) == value)]
	  vcfRaw(x) <- vcfRaw(x)[getCHROM(vcfRaw(x)) %in% value ]
	  invariantRaw(x) <-  invariantRaw(x)[seqnames(invariantRaw(x)) == value]
}
  validObject(x)
  x
})

# Extract a bed file with all sites sequenced
setGeneric("totalSites", function(x) standardGeneric("totalSites"))
#' @export
setMethod("totalSites", signature = "ParentClass", function(x) {
  total <- reduce(c(vcf2Grange(vcf(x)), invariant(x), fixed(x)))
  total_nb <- sum(width(total))
  message(sprintf("Total nb of sites sequenced : %s", total_nb))
  validObject(total)
  total
})


# ==========================================================================
#                            Initiator
# ==========================================================================

#' Create an empty filterVCF S4 object
#' 
#' This function is used to initialize an filterVCF S4 object with all the
#' necessary slots accesible with @. The slots can be listed by using
#' slotNames(filterVCF()). Initializing with this function is equivalent to
#' new("filterVCF") an ensure that the slots have the correct class.
#
#' 
#' @return empty filterVCF object
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords filterVCF
#' @examples
#' 
#' data <- filterVCF()
#' data
#' slotNames(data)
#' 
#' @export
filterVCF <- function(){
	new("filterVCF")
}

# ==========================================================================
#                             Validator
# ==========================================================================