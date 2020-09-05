# ==========================================================================
#                        Main functions
# ==========================================================================


#' Extract the META info from a vcf
#' 
#' Extract the META info from a vcfR object and can report plot.
#' 
#' 
#' @param vcf a vcfR object
#' @param vcf_file names fro plotting if opt$plot = TRUE
#' @param ... options from initialise_option()
#' @return list
#' \item{META }{META data extracted} 
#' \item{Description_META }{description taken from the vcf header} 
#' \item{INFO_per_snp }{META data per SNPs} 
#' \item{CI_Stats }{Statistic for each META data} 
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$plot <- TRUE
#' stats_META("my.vcf", vcf_file="Plot_metadata", opt)
#' 
#' @export
stats_META <- function(vcf, vcf_file="vcf", ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}

	if(nrow(vcf) == 0){
		message("\nNo site left after filtering!!")
		return(list())
	} 
	#Export some stats adapted for GATK,
	#it can also output from other caller
	#but not tested.
	META <- grep("INFO", queryMETA(vcf), value = TRUE)
	META <- gsub(".+=(.+)$", "\\1", META)
	META <- grep("AC|AF|AN", META, value = TRUE, invert =TRUE)

	Description_META <- lapply(META, function(x) paste(x, grep("Description", queryMETA(vcf, element = x)[[1]], value = TRUE), sep = ":"))
	INFO_per_snp <- lapply(META, function(x) extract.info(vcf, element = x, as.numeric = TRUE))
	names(INFO_per_snp) <- names(Description_META) <- META

	mat_INFO_per_snp <- as.data.frame(matrix(unlist(INFO_per_snp), ncol=length(INFO_per_snp)))
	index_col <-  colSums(mat_INFO_per_snp, na.rm=TRUE) != 0
	mat_INFO_per_snp <- mat_INFO_per_snp[,index_col]
	mat_INFO_per_snp <- mat_INFO_per_snp[complete.cases(mat_INFO_per_snp),]
	colnames(mat_INFO_per_snp) <- META[index_col]
	if(opt$verbose) message(sprintf("Computing stats for INFO :\n%s", paste(META[index_col], collapse = " ")), appendLF=TRUE)

	CI_Stats <- apply(mat_INFO_per_snp,2, quantile, probs = c(0.025, 0.5, 0.975))

	if(opt$plot){
		pdf(paste0(gsub(".vcf|.vcf.gz", "",vcf_file), "_", "filter_stat_overview.pdf"))
		#pdf()
		par(mfrow=c(2, 2), cex.main = 0.5)
		Hist_out <- suppressWarnings(invisible(lapply(META, function(x) try(hist(INFO_per_snp[[x]], 20, xlab = x, main=Description_META[[x]], col = "aquamarine2"),silent=TRUE))))

		pc <- try(prcomp(scale(mat_INFO_per_snp)))
		par(mfrow=c(1, 1), cex.main = 0.6)
		try(biplot(pc, choices =1:2, cex=.6, pc.biplot = TRUE))
		try(biplot(pc, choices =3:4, cex=.6, pc.biplot = TRUE))
		graphics.off()
	}
	message("Done...", appendLF=TRUE)
	return(list(META=META, Description_META=Description_META, INFO_per_snp=INFO_per_snp, CI_Stats=CI_Stats))
}


#' Save large output from a filterVCF
#' 
#' This function is used in the pipeline filterVCF.R to save output from a run
#' of filtering. The filterVCF object contains information about the sites that
#' were filtered, the invairant, the fixed, the bedfile imputed,etc. This
#' function creates a directory with the resulting filtered vcf file, invariant
#' and fixed bedfile and all the filterts that have been applied as bedfiles.
#' 
#' 
#' @param RES a filterVCF object
#' @param output_dir output directory
#' @param bed_outputdir output directory to save filters
#' @return NULL
#' 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' opt <- GetOpt(filterVCF_object)
#' save_output(filterVCF_object, output_dir = "./", bed_outputdir = "filters", opt)
#' 
#' @export
save_output <- function(RES, output_dir = "./", bed_outputdir, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}

	if(!is.na(opt$output_file)){
		vcf_file_short <- paste0(basename(gsub(".vcf$|.vcf.gz$", "", opt$output_file)))
	}else{
		vcf_file_short <- paste0(basename(gsub(".vcf$|.vcf.gz$", "", opt$vcf_file)))		
	}
	bed_window_filter_dir <- paste0(bed_outputdir, "/window_filters")
	if(!dir.exists(bed_window_filter_dir)){
		dir.create(bed_window_filter_dir)
	}

	if(!is.na(opt$bed_file[1])){
		Grange2bed(BEDFilter(RES), "Filter by bed file",   bed_window_filter_dir, paste0(vcf_file_short, "_filtered_by_", "bed_GRange_merged",  ".bed"))
		}

	if(isTRUE(opt$filter_repeats_by_windows)){
		Grange2bed(Filter_WindowsProp(RES), "windows with half repeats (or other bedfile specified in --repeats)",  bed_window_filter_dir, paste0(vcf_file_short, "_filtered_by_", "windows_with_half_repeats_to_remove",  ".bed"))
		}

	if(isTRUE(opt$filter_fix_het)){
		Grange2bed(FixHet_RefPop(RES), "Fixed contiguous heterozygous sites in 50bp windows in the population specified in --filter_fix_het_contiguous_in_pop", bed_window_filter_dir, paste0(vcf_file_short, "_filtered_by_", "fix_het_pop",  ".bed"))
		}

	if(isTRUE(!is.na(opt$filter_high_DP_standardized))){
		names(Normalized_DP_Filter(RES)[[1]]) <- names(Normalized_DP_Filter(RES)[[2]]) <- paste0("median_normDP_", parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)$threshold)

		for(x in names(Normalized_DP_Filter(RES))){
			for(y in names(Normalized_DP_Filter(RES)[[x]])){
					Grange2bed(Normalized_DP_Filter(RES)[[x]][[y]], sprintf("%s with normalized depth higher than the threshold %s", x, y), bed_window_filter_dir, paste0(vcf_file_short, "_filtered_by_", "standardized_DP_", y, ifelse(x=="total_high_DP_windows", "_per_sample", "_across_sample"),  ".bed") )
			}
	}
	}

	#Save bed file for sites removed by diffrent filters to intersect with vcf_bed_Grange
	bed_single_filter_dir <- paste0(bed_outputdir, "/individual_filters")

	filter2save <- c(
	"QD",
	"SOR",
	"MQRankSum",
	"FS",
	"MQ",
	"ReadPosRankSum",
	"InbreedingCoeff",
	"fix_het",
	"all_het",
	"bi_allelic",
	"missing",
	"indel")

	vcf_bed_GRange <- vcf2Grange(vcfRaw(RES))
	if(!dir.exists(bed_single_filter_dir)){
		dir.create(bed_single_filter_dir)
	}

	void <- lapply(filter2save, function(x) filter2bed(Filters(RES)[[x]], x, vcf_bed_GRange, bed_single_filter_dir, paste0(vcf_file_short, "_filtered_by_", x,  ".bed")))

	#save removed sites
	master_filter_Grange <- removed_sites(RES)

	if(sum(width(invariantRaw(RES))) != 0){
		master_filter_Grange <- c(master_filter_Grange, removed_sites_inv(RES))
	}

	Grange2bed(master_filter_Grange, "Removed sites", bed_outputdir, paste0(vcf_file_short, "_removed_sites",  ".bed"))

	#save fixed invariant
	if(isTRUE(length(fixed(RES)) != 0) ){
		Grange2bed(fixed(RES), "Fixed homozygous ALT sites , keep that filter if a repolarization is done later", bed_outputdir, paste0(vcf_file_short, "_ALT_fixed",  ".bed"))
	}

	#save invariant
	if(isTRUE(length(invariant(RES)) != 0)){
		Grange2bed(invariant(RES), "Invariable sites including REF/REF and missing", bed_outputdir, paste0(vcf_file_short, "_invariable",  ".bed"))
	}

	list_bed_files <- list.files(path=outputdir, pattern=".bed$", recursive = TRUE, full.names = TRUE)
	void <- lapply(list_bed_files, gzip, overwrite=TRUE)
	#gzip all bed file
	#write vcf file
	if(opt$verbose) message(sprintf("Save filtered vcf to %s: ", paste0(vcf_file_short, "_filtered",  ".vcf.gz")))
	write.vcf(vcf(RES), file =paste0(vcf_file_short, "_filtered",  ".vcf.gz"))
}



#' Main wrapper to filter a VCF file
#' 
#' This function is the main core of the pipeline filterVCF.R. It will perform
#' a series of filtering steps depending of the options provided in a opt list
#' that can be initiated with initialise_option(). This function is not expected
#' to be used directly by the user which should first run the pipeline
#' filterVCF.R on a complete non filtered VCF (including invariant). 
#' 
#' 
#' see the help of filterVCF.R for more detail
#' Rscript filterVCF.R --help
#' 
#' 
#' @param vcf_file path to vcf file
#' @param ... options from initialise_option(). This parameter regulates all the filtering steps that will be applied to the data.
#' @return a filterVCF object.
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~main
#' @examples
#' 
#' opt <- parse_args(OptionParser(option_list=option_list))
#' main_filter_VCF("my.vcf",  opt)
#' 
#' @export
main_filter_VCF <- function(vcf_file, ...){
	.libPaths("/proj/uppstore2018024/private/Rpackages/")
	suppressPackageStartupMessages(suppressMessages(try(library(filterVCF))))
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}

	vcf <- read.vcfR(vcf_file, limit = 1e+08, verbose = opt$verbose, convertNA = FALSE)
	vcf@fix[,"ID"] <- NA
	
	if(!is.na(opt$keep_sites)){
		sites2keep <- bed2Grange(opt$keep_sites)
		vcf_bed_GRange <- vcf2Grange(vcf)
		index_sites2keep <- findOverlaps(vcf_bed_GRange, sites2keep, minoverlap=1, ignore.strand=TRUE)
		vcf <- vcf[from(index_sites2keep)]
		if(opt$verbose) message(sprintf("Keeping only sites interesecting with %s\nResults in %s sites before filtering", opt$keep_sites, nrow(vcf)), appendLF=TRUE)
	}
	#select only some chromosome
	#return a empty RES if the chr i not present in the vcf
	if(!is.na(opt$chr[1])){
		vcf <- vcf[getCHROM(vcf) %in% opt$chr]
		if(nrow(vcf) == 0 ){
			return(as.filterVCF(list(
						vcf=vcf,
						vcf_filtered=vcf,
						vcf_inv=GRanges(),
						vcf_inv_filtered=GRanges(),
						fixed_inv_Grange=GRanges(),
						reference=new("DNAStringSet"),
						master_filter=logical(0),
						master_filter_inv=logical(0),
						filters=list(logical(0)),
						filters_inv=list(logical(0)),
						removed_sites=GRanges(),
						removed_sites_inv=GRanges(),
						bed_GRange_merged=GRanges(),
						windows_with_half_repeats_to_remove=GRanges(),
						fix_het_pop=GRanges(),
						filter_DP_repeats=list(high_DP_windows_means_across_sample=GRanges(), total_high_DP_windows=GRanges()),
						stats_vcf_filtered=NA_character_,
						opt = opt)))
		}
	}

	if(!is.na(opt$sample)){
		opt$sample_list <- parse_opt_sample(opt$sample)
		vcf <- vcf[,c("FORMAT", opt$sample_list)]
		warnings("\nInformation in the INFO field will not be recalculated, use GATK -T SelectVariants on the output vcf to recreate them.\n")
	}
	# extract_genotype
	gt <- extract_gt_ff(vcf)

	#remove sample with too much missing data
	if(!is.na(opt$missing_per_individual)){
		passing_missing_ind <- missing_per_individual(gt, opt$missing_per_individual)
		vcf <- vcf[,c("FORMAT", passing_missing_ind)]
		gt <- gt[passing_missing_ind]
	}

	#Start by filtering the invariant out and keep the vcf for further saving
	filter_invariant <- create_invariant_filter(vcf, gt=gt, opt)
	gc(verbose=FALSE)
	vcf_inv <- vcf[filter_invariant]
	vcf <- vcf[!filter_invariant]

	if(opt$verbose) message("Done...", appendLF=TRUE)
	#explore the filtering option in the metadata and plot some reports. If the vcf is too big start by subsampling 100000
	if(!is.na(opt$split) & isTRUE(opt$split != 1)){
		#dont calculate stats on splitted file
		#it will be done later
		# note that is a subsample is taken the info will not match for the sample
		# it is better to use SelectVariants first
		if(nrow(vcf) > 1e5){
			if(opt$verbose) message(sprintf("The VCF has %s variants, taking a subsample of 100000 to draw stats plot", nrow(vcf)), appendLF=TRUE)
			vcf_for_stats <- vcf[sample(1:nrow(vcf), 1e5),]
		}else{
			vcf_for_stats <- vcf
		}

		#get some stats
		stats_vcf <- suppressWarnings(stats_META(vcf_for_stats, vcf_file=vcf_file, opt))
	}

	# extract_genotype
	gt <- gt[as.ff(!filter_invariant),]
	#filter DP per individual, if used filter missing should be recalculated
	#as we will replace the genotype by NA is DP is < threshold > .
	# vcf_bak <- vcf
	# vcf <- vcf_bak
	if(!is.na(opt$filter_depth)){
		opt$filter_depth_parsed <- parse_filter_depth(opt$filter_depth)
		vcf <- change_gt_by_DP_filter(vcf, min_DP=opt$filter_depth_parsed$min_DP, max_DP=opt$filter_depth_parsed$max_DP, opt)
		#We need to re extract the genotype
		rm(gt)
		gt <- extract_gt_ff(vcf)

		#After changing the genotyoe to NA some sites can become invariant
		#create a new  filter for invariant and
		#adding the new inv to the old while filtering the snp
		#it works even if no sites is filtered
		filter_invariant_dp <- create_invariant_filter(vcf, gt=gt, opt)

		vcf_inv <- rbind2(vcf_inv,  vcf[filter_invariant_dp])
		vcf <- vcf[!filter_invariant_dp]
		gt <- gt[as.ff(!filter_invariant_dp),]
		#re add the invariant
	}

	if(!is.na(opt$filter_genotype_quality)){
		opt$filter_genotype_quality <- parse_filter_GQ(opt$filter_genotype_quality)
		
		vcf <- change_gt_by_GQ_filter(vcf, invariable = FALSE, GQ_threshold=opt$filter_genotype_quality$GQ, opt)
		#We need to re extract the genotype
		rm(gt)
		gt <- extract_gt_ff(vcf)

		#After changing the genotyoe to NA some sites can become invariant
		#create a new  filter for invariant and
		#adding the new inv to the old while filtering the snp
		#it works even if no sites is filtered
		filter_invariant_dp <- create_invariant_filter(vcf, gt=gt, opt)

		vcf_inv <- rbind2(vcf_inv,  vcf[filter_invariant_dp])
		vcf <- vcf[!filter_invariant_dp]
		gt <- gt[as.ff(!filter_invariant_dp),]
		#re add the invariant
	}

	#save the inv to temp file to retrieve after to sve memory usage
	tmp_vcf_inv <- tempfile()
	saveRDS(vcf_inv, file = tmp_vcf_inv)
	rm(vcf_inv)

	filters <- list()
	# Filter the variant
	# GATK recommendation modified after exploring the data. It is a bit more stringent
	# and try to remove excess of hets due to wrong call
	# "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	# "QD < 5.0 || FS > 60.0 || MQ < 50.0 || MQRankSum < -5 || ReadPosRankSum < -4.0 || SOR > 3 || InbreedingCoeff < -0.2" \
	if(!is.null(opt$filter_GATK_info)){
		filter_GATK_info <- parse_GATK_filter_option(opt$filter_GATK_info)

		for(i in names(filter_GATK_info)){
			filters[[i]] <- create_filter_META(vcf, i, filter_GATK_info[[i]])
		}
	}

	# filter for fix het
	if(opt$filter_fix_het){
		filters[["fix_het"]] <- create_filter_fix_het(gt, opt)
	}

	# filter for fix het
	if(opt$filter_all_het){
		filters[["all_het"]] <- create_filter_all_het(gt, opt)
	}

	#filter for bi allelic
	if(opt$biallelic){
		if(opt$verbose) message("Apply filter to keep bi-allelic sites", appendLF=TRUE)
		filters[["bi_allelic"]] <- !is.biallelic(vcf)
	}

	#filter for missing
	if(!is.na(opt$missing)){
		filters[["missing"]] <- create_filter_missing(gt, opt$missing)
	}

	# filter snps only remove indels
	if(opt$filter_indel){
		filters[["indel"]] <- create_indel_filter(vcf, opt)
	}

	# Reading all bedfiles into Grange object
	if(!is.na(opt$bed_file[1])){
		if(opt$verbose) message("Apply filter from bedfiles", appendLF=TRUE)

		bed_GRange <- lapply(opt$bed_file, bed2Grange)

		# Merging Grange object into one to intersect with the vcf
		bed_GRange_merged <- suppressWarnings(do.call("c",  bed_GRange))

		# Converting vcf to Grange
		vcf_bed_GRange <- vcf2Grange(vcf)

		#finding the overlap between the vcf and the bed file (bedtools subtract)
		filters[["filter_bed"]] <- !is.na(GenomicRanges::findOverlaps(vcf_bed_GRange, bed_GRange_merged , select ="first", ignore.strand=TRUE))
	}else{
		bed_GRange_merged <- NULL
	}

	if(isTRUE(opt$filter_repeats_by_windows)){
		#repeats windows with more than 50%
		windows_with_half_repeats_to_remove <- suppressWarnings(create_filter_repeats_in_windows(opt$reference, opt$repeats, vcf, 20000, 0.5))
		filters[["filter_repeats_by_windows"]] <- !is.na(GenomicRanges::findOverlaps(vcf_bed_GRange, windows_with_half_repeats_to_remove , select ="first", ignore.strand=TRUE))
	}else{
		windows_with_half_repeats_to_remove <- NULL
	}

	if(isTRUE(!is.na(opt$filter_high_DP_standardized))){
		#custom filter on DP per windows
		opt$filter_high_DP_standardized <- parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)

		filter_DP_repeats <- create_filter_high_DP_standardized(opt$reference, vcf, threshold = opt$filter_high_DP_standardized$threshold, threshold_CI = opt$filter_high_DP_standardized$threshold_CI, windows_size= opt$filter_high_DP_standardized$windows_size, overlapping= opt$filter_high_DP_standardized$slidding, percent_passing_filter = opt$filter_high_DP_standardized$percent_passing_filter, opt)
		gc()
		filters[["filter_DP_repeats"]] <- !is.na(GenomicRanges::findOverlaps(vcf_bed_GRange, filter_DP_repeats[[2]][[1]] , select ="first", ignore.strand=TRUE))
	}else{
		filter_DP_repeats <- NULL
	}

	if(!is.na(opt$filter_fix_het_contiguous_in_pop)){
		#filter cluster of fix het in the swedish pop
		fix_het_pop <- create_filter_fix_het_contiguous_in_pop(gt, vcf, opt$filter_fix_het_contiguous_in_pop, windows_size = 50, max_nb_contigous_hets = 1)
		vcf_bed_GRange <- vcf2Grange(vcf)
		filters[["filter_pop_het_contiguous"]] <- !is.na(GenomicRanges::findOverlaps(vcf_bed_GRange, fix_het_pop , select ="first", ignore.strand=TRUE))
	}else{
		fix_het_pop <- NULL
	}

	# Sum up all filters to create a master filter
	master_filter <- Reduce("+", filters) > 0
	if(length(master_filter) == 0) master_filter <- rep(FALSE, nrow(vcf))
	# #information on the number of snps remove by each filters
	# snp_remove_by_filter <- lapply(filters, sum, na.rm = TRUE)

	#write a report of the snps that were filtered

	#Remove the fixed sites in the vcf and it to the invariant vcf.
	#Save a bedfile with the fixed sites
	filters[["inv_fix_sites"]] <- create_filter_fixed_sites(gt)
	filters[["inv_fix_sites"]] <- !master_filter & filters[["inv_fix_sites"]]
	fixed_inv_Grange <- vcf2Grange(vcf[filters[["inv_fix_sites"]],])

	rm(gt)
	gc()

	fix_vcf <- vcf[filters[["inv_fix_sites"]],]
	#join the invariant with the fixed in a vcf
	vcf_inv <- rbind2(readRDS(tmp_vcf_inv), fix_vcf)

	#filter the vcf file with all the filter at once.
	#single filter will be outputed to a folder
	#in order to create removed sites bedfile or
	#refiltered the original vcf using only some
	#filters with bedtools subtract.
	vcf_filtered <- vcf[!master_filter & !filters[["inv_fix_sites"]],]
	stats_vcf_filtered <- stats_META(vcf_filtered, vcf_file=paste0(vcf_file, "_filtered"), opt)

	# ==========================================================================
	#       apply some filters to vcf_inv
	# ==========================================================================
	
	if(!is.na(opt$filter_depth)){
		vcf_inv <- change_gt_by_DP_filter(vcf_inv, min_DP=opt$filter_depth_parsed$min_DP, max_DP=opt$filter_depth_parsed$max_DP, opt)
	}
	
	#filter on Ref Genotype Quality
	if(!is.na(opt$filter_genotype_quality)){
		vcf_inv <- change_gt_by_GQ_filter(vcf_inv, invariable = TRUE, GQ_threshold=opt$filter_genotype_quality$RGQ, opt)
	}

	filters_inv <- list()

	#filter for bi allelic
	if(opt$biallelic){
		filters_inv[["bi_allelic"]] <- !is.biallelic(vcf_inv)
	}

	#filter for missing
	if(!is.na(opt$missing)){
		filters_inv[["missing"]] <- create_filter_missing(vcf_inv, opt$missing)
	}

	# filter snps only remove indels
	if(opt$filter_indel){
		filters_inv[["indel"]] <- create_indel_filter(vcf_inv, invariable=TRUE)
	}

	# Converting vcf to Grange
	vcf_inv_bed_GRange <- vcf2Grange(vcf_inv)

	#finding the overlap between the vcf and the bed file (bedtools subtract)
	if(!is.na(opt$bed_file[1])){
		#finding the overlap between the vcf and the bed file (bedtools subtract)
		filters_inv[["filter_bed"]] <- !is.na(GenomicRanges::findOverlaps(vcf_inv_bed_GRange, bed_GRange_merged , select ="first", ignore.strand=TRUE))
	}

	#repeats windows with more than 50%
	if(opt$filter_repeats_by_windows){
		filters_inv[["filter_repeats_by_windows"]] <- !is.na(GenomicRanges::findOverlaps(vcf_inv_bed_GRange, windows_with_half_repeats_to_remove , select ="first", ignore.strand=TRUE))
	}

	if(isTRUE(!is.na(opt$filter_high_DP_standardized))){
		#custom filter on DP per windows
		vcf_inv_bed_GRange <- vcf2Grange(vcf_inv)
		filter_DP_repeats_inv <- create_filter_high_DP_standardized(opt$reference, vcf_inv, threshold = opt$filter_high_DP_standardized$threshold, threshold_CI = opt$filter_high_DP_standardized$threshold_CI, windows_size= opt$filter_high_DP_standardized$windows_size, overlapping= opt$filter_high_DP_standardized$slidding, opt)
		gc()
		filters_inv[["filter_DP_repeats"]] <- !is.na(GenomicRanges::findOverlaps(vcf_inv_bed_GRange, filter_DP_repeats_inv[[1]][[1]] , select ="first", ignore.strand=TRUE))
	}else{
		filter_DP_repeats_inv <- NULL
	}


	# Sum up all filters to create a master filter
	master_filter_inv <- Reduce("+", filters_inv) > 0

	#information on the number of snps remove by each filters
	#inv_remove_by_filter <- lapply(filters_inv, sum, na.rm = TRUE)

	vcf_inv_filtered <- vcf2Grange(vcf_inv[!master_filter_inv], metadata = "present")
	###############
	vcf_bed_GRange <- vcf2Grange(vcf)
	removed_sites <- reduce(vcf_bed_GRange[master_filter])
	removed_sites_inv <- reduce(vcf_inv_bed_GRange[master_filter_inv])

	RES <- as.filterVCF(list(
			vcf=vcf,
			vcf_filtered=vcf_filtered,
			vcf_inv=vcf2Grange(vcf_inv),
			vcf_inv_filtered=vcf_inv_filtered,
			fixed_inv_Grange=fixed_inv_Grange,
			reference=readDNAStringSet(opt$reference),
			master_filter=master_filter,
			master_filter_inv=master_filter_inv,
			filters=filters,
			filters_inv=filters_inv,
			removed_sites=removed_sites,
			removed_sites_inv=removed_sites_inv,
			bed_GRange_merged=bed_GRange_merged,
			windows_with_half_repeats_to_remove=windows_with_half_repeats_to_remove,
			fix_het_pop=fix_het_pop,
			filter_DP_repeats=filter_DP_repeats,
			stats_vcf_filtered=NA,
			opt=opt))
			#stats_vcf_filtered=stats_vcf_filtered$CI_Stats))

	rm(vcf_inv)
	rm(vcf)
	gc()

	if(isTRUE(opt$debug)){
		if(!is.na(opt$output_file)){
			saveRDS(RES, paste0(getwd(), "/", basename(paste0(gsub(".vcf$|.vcf.gz$", "", opt$output_file), ".rds"))))		
		}else{
			saveRDS(RES, paste0(getwd(), "/", basename(gsub(".vcf$|.vcf.gz$", ".rds", vcf_file))))			
		}
	}
	if(opt$verbose) Sys.procmem()

	#remove vcf_inv and inv_filtered to save some disk space?
	return(RES)
}

