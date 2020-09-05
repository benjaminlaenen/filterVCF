
# ==========================================================================
#                SET Class VIFobject
# ==========================================================================
#' An S4 class to for VariantInvariantFixed (VIF).
#'
#' @slot vcf_filtered The filtered vcf of class "vcfR",
#' @slot vcf_inv_filtered The list of invariant filtered sites of class "GRanges",
#' @slot fixed_inv_Grange The list of invariant fixed sites of class"GRanges",
#' @slot reference A DNAStringSet object with the reference used to create the object,
#' @slot removed_sites All sites removed in snp of class "GRanges",
#' @slot removed_sites_inv All sites removed in invariant of class "GRanges",
#' @slot Divergence from outgroupd of class "GRanges, structuee is chr pos REF out1 out2 ... outn
#' @slot Repolarization key of class "GRanges", produce using either likelohood (Keightley) or parsimony. 
#' @slot bootstrap List to store bootstrap index to reproduce boostrapped dataset quickly,
#' @slot opt List of option used by the pipeline for VIF,
#' @export
#' @import IRanges GenomicRanges rtracklayer Biostrings
setClass("VIF", slots = c(
vcf_filtered                        = "vcfR",
vcf_inv_filtered                    = "GRanges",
fixed_inv_Grange                    = "GRanges",
reference                           = "DNAStringSet",
removed_sites                       = "GRanges",
removed_sites_inv                   = "GRanges",
divergence                          = "GRanges",
alternate                           = "GRanges",
repolKey                            = "GRanges",
bootstrap                           = "list",
opt                                 = "list"),
contains = "ParentClass"
)

setGeneric("show.VIF", function(object) {
  standardGeneric("show.VIF")
})

#' An S4 methods for object filterVCF
#'
#' Print information on the current object.
#' @export
setMethod("show", signature(object = "VIF"), function(object) {
  nsample <- ncol(vcf(object)@gt) - 1
  n_SNP <- nrow(vcf(object))
  SNP <- n_SNP
  n_inv <- sum(width(invariant(object)))
  n_fixed <- sum(width(fixed(object)))
  divergence <- sum(width(divergence(object)))
  printTable <- c(SNPs=SNP, Invariable =n_inv, Fixed = n_fixed)
  cat("
  # ==========================================================================
  #                         Object of class VIF
  # ==========================================================================\n")
  cat(paste("The number of sample is : ", nsample))
  cat("\n\n")
  cat("Number of sites")
  cat("\n\n")
  print(printTable)
  cat("\n")
  cat("
  # ==========================================================================
  #                        
  # ==========================================================================\n")
})


# ==========================================================================
#                       Generic
# ==========================================================================
setGeneric("as.VIF", function(object) standardGeneric("as.VIF"))
setGeneric("as.list.VIF", function(object) standardGeneric("as.list.VIF"))
setGeneric("repairVIF", function(object) standardGeneric("repairVIF"))

# ==========================================================================
#                         Methods
# ==========================================================================

#' An S4 methods for object list
#'
#' Coerce a list to VIF object. Some slots can be empty.
#' @export
setMethod("as.VIF", signature(object = "list"), function(object) {
  res <- VIF()
  res@vcf_filtered                               = object$vcf_filtered
  res@vcf_inv_filtered                           = object$vcf_inv_filtered
  res@fixed_inv_Grange                           = object$fixed_inv_Grange
  res@removed_sites                              = object$removed_sites
  res@removed_sites_inv                          = object$removed_sites_inv
  res@divergence                                 = if(is.null(object$divergence)) GRanges() else object$divergence
  res@repolKey                                   = if(is.null(object$repolKey))GRanges() else object$repolKey
  res@alternate                                  = if(is.null(object$alternate))GRanges() else object$alternate
  res@reference                                  = if(is.null(object$reference)) new("DNAStringSet") else object$reference
  res@bootstrap                                  = if(is.null(object$bootstrap)) initialise_bootstrap_list() else object$bootstrap
  res@opt                                        = if(is.null(object$opt)) parse_args(OptionParser(option_list=initialise_option())) else object$opt
  return(res)
})


#' An S4 methods for object filterVCF
#'
#' Coerce a filterVCF to VIF object. Some slots can be empty.
#' @export
setMethod("as.VIF", signature(object = "filterVCF"), function(object) {
  res <- VIF()
  res@vcf_filtered                               = object@vcf_filtered
  res@vcf_inv_filtered                           = object@vcf_inv_filtered
  res@fixed_inv_Grange                           = object@fixed_inv_Grange
  res@removed_sites                              = object@removed_sites
  res@removed_sites_inv                          = object@removed_sites_inv
  res@divergence                                 = GRanges()
  res@repolKey                                   = GRanges()
  res@alternate                                  = get_alternate_GRange(object)
  res@reference                                  = object@reference
  res@bootstrap                                  = initialise_bootstrap_list()
  res@opt                                        = if(!.hasSlot(object, "opt")) parse_args(OptionParser(option_list=initialise_option())) else object@opt
  return(res)
})



#' An S4 methods for object VIF
#'
#' Coerce a VIF object to list object. If slots are empty they are filled up by object of the right class.
#' @export
setMethod("as.list.VIF", signature(object = "VIF"), function(object) {
  res <- list()
  res$vcf_filtered                               = if(.hasSlot(object, "vcf_filtered")) object@vcf_filtered else new("vcfR") 
  res$vcf_inv_filtered                           = if(.hasSlot(object, "vcf_inv_filtered")) object@vcf_inv_filtered else new("GRanges")
  res$fixed_inv_Grange                           = if(.hasSlot(object, "fixed_inv_Grange")) object@fixed_inv_Grange else new("GRanges")
  res$reference                                  = if(.hasSlot(object, "reference")) object@reference else new("DNAStringSet")
  res$removed_sites                              = if(.hasSlot(object, "removed_sites")) object@removed_sites else new("GRanges")
  res$removed_sites_inv                          = if(.hasSlot(object, "removed_sites_inv")) object@removed_sites_inv else new("GRanges")
  res$divergence                                 = if(.hasSlot(object, "divergence")) object@divergence else new("GRanges")
  res$repolKey                                   = if(.hasSlot(object, "repolKey")) object@repolKey  else new("GRanges")
  res$alternate                                  = if(.hasSlot(object, "alternate")) object@alternate  else new("GRanges")
  res$opt                                        = if(.hasSlot(object, "opt")) object@opt else parse_args(OptionParser(option_list=initialise_option()))
  res$bootstrap                                  = if(.hasSlot(object, "bootstrap")) object@opt else initialise_bootstrap_list()
  res$stats_vcf_filtered                         = if(.hasSlot(object, "stats_vcf_filtered")) object@stats_vcf_filtered else NA_character_
  return(res)
})


#' An S4 methods for object VIF
#'
#' Repair an object by filling the missing slot with empty object. Function for compatibility if the definition changes between version.
#' @export
setMethod("repairVIF", signature(object = "VIF"), function(object) {
  return(as.VIF(as.list(object)))
  }
)


#' An S4 methods for object list
#'
#' Initialise an empty VIF object
#' @export
#' 
setMethod("initialize", "VIF", 
    function(.Object, ...) {
      .Object <- callNextMethod()
      .Object@opt <- parse_args(OptionParser(option_list=initialise_option()), args = NULL)
      .Object@bootstrap <- initialise_bootstrap_list()
      .Object
    }
)


#' An S4 methods for object VIF
#'
#' Coerce a list to filterVCF object. Some slots can be empty.
#' @export
setMethod("as.filterVCF", signature(object = "VIF"), function(object) {
  res <- as.filterVCF(as.list(object))
  return(res)
  }
)


# ==========================================================================
#                            Initiator
# ==========================================================================


#' Create a empty VIF object
#' 
#' Create a empty VIF object
#' 
#' 
#' @return VIF object
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  put references to the literature/web site here
#' @keywords ~initiator
#' @examples
#' 
#' VIF()
#' 
#' @export
VIF <- function(...){
  new("VIF", ...)
}

# ==========================================================================
#                          Accessor
# ==========================================================================

# Acces the divergence
setGeneric("divergence", function(x) standardGeneric("divergence"))
#' @export
setMethod("divergence", signature = "VIF", function(x) x@divergence)

setGeneric("divergence<-", function(x, value) standardGeneric("divergence<-"))
#' @export
setMethod("divergence<-", signature = "VIF", function(x, value) {
  x@divergence <- value
  validObject(x)
  x
})

# Acces the alternate
setGeneric("alternate", function(x) standardGeneric("alternate"))
#' @export
setMethod("alternate", signature = "VIF", function(x) x@alternate)

setGeneric("alternate<-", function(x, value) standardGeneric("alternate<-"))
#' @export
setMethod("alternate<-", signature = "VIF", function(x, value) {
  x@alternate <- value
  validObject(x)
  x
})

# Acces the repolarization key
setGeneric("repolKey", function(x) standardGeneric("repolKey"))
#' @export
setMethod("repolKey", signature = "VIF", function(x) x@repolKey)

setGeneric("repolKey<-", function(x, value) standardGeneric("repolKey<-"))
#' @export
setMethod("repolKey<-", signature = "VIF", function(x, value) {
  x@repolKey <- value
  validObject(x)
  x
})


# Acces the bootstrap matrix
setGeneric("bootstrapMat", function(x) standardGeneric("bootstrapMat"))
#' @export
setMethod("bootstrapMat", signature = "VIF", function(x) x@bootstrap$bootstrap_matrix)

setGeneric("bootstrapMat<-", function(x, value) standardGeneric("bootstrapMat<-"))
#' @export
setMethod("bootstrapMat<-", signature = "VIF", function(x, value) {
  x@bootstrap$bootstrap_matrix <- value
  validObject(x)
  x
})

# Acces and set the number of bootstrap
setGeneric("Nboot", function(x) standardGeneric("Nboot"))
#' @export
setMethod("Nboot", signature = "VIF", function(x) x@bootstrap$nb_bootstrap)

setGeneric("Nboot<-", function(x, value) standardGeneric("Nboot<-"))
#' @export
setMethod("Nboot<-", signature = "VIF", function(x, value) {
  x@bootstrap$nb_bootstrap <- value
  validObject(x)
  x
})

# Acces and set the how bootstrap is performed
setGeneric("bootBy", function(x) standardGeneric("bootBy"))
#' @export
setMethod("bootBy", signature = "VIF", function(x) x@bootstrap$bootstrap_by)

setGeneric("bootBy<-", function(x, value) standardGeneric("bootBy<-"))
#' @export
setMethod("bootBy<-", signature = "VIF", function(x, value) {
  x@bootstrap$bootstrap_by <- value
  validObject(x)
  x
})

# Acces and set the seed(s) for bootstrap 
setGeneric("bootSeed", function(x) standardGeneric("bootSeed"))
#' @export
setMethod("bootSeed", signature = "VIF", function(x) x@bootstrap$seed)

setGeneric("bootSeed<-", function(x, value) standardGeneric("bootSeed<-"))
#' @export
setMethod("bootSeed<-", signature = "VIF", function(x, value) {
  x@bootstrap$seed <- value
  validObject(x)
  x
})

# Acces the Granges defining the bootstrap interval
# either custom range (gff, gene length) or windows 
setGeneric("bootGRange", function(x) standardGeneric("bootGRange"))
#' @export
setMethod("bootGRange", signature = "VIF", function(x) x@bootstrap$intersect_Granges)

setGeneric("bootGRange<-", function(x, value) standardGeneric("bootGRange<-"))
#' @export
setMethod("bootGRange<-", signature = "VIF", function(x, value) {
  x@bootstrap$intersect_Granges <- value
  validObject(x)
  x
})

