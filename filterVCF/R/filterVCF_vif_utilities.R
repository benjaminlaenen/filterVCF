
## ==========================================================================
#         Utility function for VIF object and produce output
# ==========================================================================
#' @export
initialise_bootstrap_list <- function(){
	list(intersect_Granges=GRanges(), bootstrap_matrix=matrix(), all_ranges_vcf= list(), all_ranges_inv = list(), all_ranges_fixed = list(), nb_bootstrap = NULL, bootstrap_by = "gene", seed = NULL)
}



#' Subsample a VIF object
#' 
#' Extract a subsample using a list of individuals and recalculate SNps, invariant and fixed sites
#' 
#' 
#' @inheritParams repolarize_VIF
#' @param subsample regex pattern to grep the VCF column names. Can be a list of individual separated by "|"
#' @param missing_allowed Percentage of missing data allowed in the subsample
#' @param gt a genotype matrix extracted using extract_gt_ff(), speeds up the function
#' @return a subsampled VIF onject
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~VIF
#' @examples
#' 
#' sub_vif  <- subsample_VIF(vif_object, subsample = "ind1|ind3|ind234")
#' 
#' @export
subsample_VIF <- function(vif, subsample, missing_allowed = 0.2, gt = NULL){
	#subsample = "ITM2"
	all_sample <- colnames(vcf(vif)@gt)[-1]
	sample_list <- sort(grep(subsample, all_sample, value = TRUE))
	VCF <- vcf(vif)[,c("FORMAT", sample_list)]
	if(is.null(gt)){
		gt <- extract_gt_ff(VCF)
	}else{
		gt <- gt[,sample_list]
	}

	#remove site that have become all missin
	filter_missing <- create_filter_missing(gt, allowed_missing_threshold =  if(is.na(missing_allowed)) GetOpt(vif)$missing else missing_allowed)

	filter_invariant <- create_invariant_filter(VCF, gt=gt)
	vcf_inv_Grange <- vcf2Grange(VCF[filter_invariant & !filter_missing])

	filter_fix_sites <- create_filter_fixed_sites(gt)
	fixed_inv_Grange <- vcf2Grange(VCF[filter_fix_sites & !filter_missing,])

	#create new result object
	res_vif <- VIF()

	invariant(res_vif) <- c(invariant(vif), vcf_inv_Grange)
	fixed(res_vif) <- c(fixed(vif), fixed_inv_Grange)

	if(any(unique(width(fixed(res_vif))) != 1)){
		fixed(res_vif) <- expand_all_pos_Granges(fixed(res_vif))
	}
	
	#remove the fixed from the invariant so the are exclusive
	invariant(res_vif) <- intersect_Granges(invariant(res_vif), fixed(res_vif), invert = TRUE)

	removed_sites <- vcf2Grange(VCF[!filter_invariant & !filter_fix_sites & filter_missing,])
	removed_sites(res_vif) <- c(removed_sites(vif), removed_sites)

	removed_sites_inv <- vcf2Grange(VCF[filter_invariant | filter_fix_sites & filter_missing,])
	removed_sites_inv(res_vif) <- c(removed_sites_inv(vif), removed_sites_inv)

	vcf(res_vif) <- VCF[!filter_invariant & !filter_fix_sites & !filter_missing]

	reference(res_vif) <- reference(vif)
	divergence(res_vif) <- intersect_Granges(divergence(vif), c(vcf2Grange(vcf(res_vif)), invariant(res_vif), fixed(res_vif) ), merge_interval = FALSE)
	alternate(res_vif) <- intersect_Granges(alternate(vif), c(vcf2Grange(vcf(res_vif)), invariant(res_vif), fixed(res_vif) ), merge_interval = FALSE)
	GetOpt(res_vif) <- GetOpt(vif)
	repolKey(res_vif) <- intersect_Granges(repolKey(vif), c(vcf2Grange(vcf(res_vif)), invariant(res_vif), fixed(res_vif)), merge_interval = FALSE)

	message(sprintf("\nKeeping %s samples:\n%s", length(sample_list), paste(sample_list, collapse = "\n")))
	message(sprintf("\nRemoving %s sites with missing data in the subsample", sum(filter_missing)))
	message(sprintf("\n%s sites became invariant in the subsample", sum(filter_invariant)))
	message(sprintf("\n%s sites became fixed in the subsample", sum(filter_fix_sites)))
	message(sprintf("\nRemaining %s variant sites in the subsampled vcf", nrow(vcf(res_vif))))
	return(res_vif)
}

#' Intersect a VIF with a GRange
#' 
#' Simple intersection. All the slots in the VIF object will be intersect with the GRange.
#' 
#' @inheritParams repolarize_VIF
#' @param Grange GRange object to intersect with
#' @param check_overlap check for overlapping ranges in the Granges. Such overlapped ranges can created duplicated entries in the VIF.
#' @param invert subtract the Grange from the VIF
#' @return a VIF object
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~VIF
#' @examples
#' 
#' 
#' 
#' @export
intersect_VIF_Grange <- function(vif, Grange, check_overlap = TRUE, invert = FALSE){
	if(class(Grange) == "character"){
		Grange <- bed2Grange(Grange)
	}

	#test if some range are overlapping
	#which will result in having duplicated sites
	#in the vcf
	if(check_overlap){
		Grange <- remove_overlapping_range(Grange)
	}

	res <- VIF()
	vcf(res) <- intersect_Granges(vcf(vif), Grange, merge_interval = FALSE, invert = invert)
	invariant(res) <- intersect_Granges(invariant(vif), Grange, merge_interval = FALSE, invert = invert)
	fixed(res) <- intersect_Granges(fixed(vif), Grange, merge_interval = FALSE, invert = invert)
	reference(res) <- reference(vif)
	alternate(res) <- intersect_Granges(alternate(vif), Grange, merge_interval = FALSE, invert = invert)
	removed_sites(res) <- intersect_Granges(removed_sites(vif), Grange, merge_interval = FALSE, invert = invert)
	removed_sites_inv(res) <- intersect_Granges(removed_sites_inv(vif), Grange, merge_interval = FALSE, invert = invert)
	divergence(res) <- intersect_Granges(divergence(vif), Grange, merge_interval = FALSE, invert = invert)
	repolKey(res) <- intersect_Granges(repolKey(vif), Grange, merge_interval = FALSE, invert = invert)
	GetOpt(res) <- GetOpt(vif)
	#keep track of the sites with used for the intersection
	GetOpt(res)$intersect_Granges <- Grange
	bootGRange(res) <- Grange
	return(res)
}

#' Split a VIF object into ranges
#' 
#' This function will split all slots in the VIF object into the number of ranges provided by the Granges.
#' 
#' 
#' @inheritParams intersect_VIF_Grange
#' @return a list of VIF object the length of Grange
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  put references to the literature/web site here
#' @examples
#' 
#' length(Grange)
#' splitted_vif  <- split_VIF(vif_object, Grange)
#' identical(length(splitted_vif), length(Grange))
#' 
#' @export
split_VIF <- function(vif, Grange){
	elementMetadata(Grange) <- NULL
	vcf_split <- split_vcfR(vcf(vif), Grange)
	invariant_split <- split_Grange(invariant(vif), Grange)
	fixed_split <- split_Grange(fixed(vif), Grange)

	removed_sites_split <- split_Grange(removed_sites(vif), Grange)
	removed_sites_inv_split <- split_Grange(removed_sites_inv(vif), Grange)
	divergence_split <- split_Grange(divergence(vif), Grange)
	alternate_split <- split_Grange(alternate(vif), Grange)
	repolKey_split <- split_Grange(repolKey(vif), Grange)
	
	vif_split <- lapply(seq_along(vcf_split), function(i) VIF(vcf_filtered = vcf_split[[i]], 
															  vcf_inv_filtered = invariant_split[[i]], 
															  fixed_inv_Grange = fixed_split[[i]],
															  removed_sites = removed_sites_split[[i]],
															  removed_sites_inv = removed_sites_inv_split[[i]],
															  divergence = divergence_split[[i]],
															  alternate = alternate_split[[i]],
															  repolKey = repolKey_split[[i]],
															  reference = reference(vif),
															  bootstrap = list(),
															  opt = GetOpt(vif)))
}
