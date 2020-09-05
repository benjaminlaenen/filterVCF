
# ==========================================================================
#              Popstat VCF
# ==========================================================================
#' Downsample site frequency spectrum
#' 
#' This function uses the formula from Nielsen (2004) to downsample a site
#' frequency spectrum (SFS) to any number of chromosome. This can be applied to
#' unfolded ot folded spectra.
#' 
#' Note the function return a rounded SFS, hence the number of sites differ sensibly before and after downsampling.
#' 
#' @param SFS an integer vector 
#' @param target_sample the final number of chromosomes desired
#' @param include_invariable Does the SFS contains invariable, 0-bin
#' @return A downsampled SFS of length of the target_sample.
#' @author  ~~Benjamin Laenen~~
#'
#' @references  Nielsen, Rasmus. "Population genetic analysis of ascertained SNP data." Human genomics 1.3 (2004): 218.
#' @keywords ~utilities
#' @examples
#' 
#' SFS <- getSFS_VIF(VIF)$SFS 
#' downsample_SFS(SFS, 10)
#' 
#' @export
downsample_SFS <- function(SFS, target_sample, include_invariable = TRUE){
	#methods from Nielsen et al 2005 for downsampling
	if(!include_invariable) SFS <- c(0, SFS)
	nj = length(SFS) - 1
	if(nj == target_sample){
		return(SFS)
	}else{
		x_down <- matrix(NA, target_sample+1, nj+1)
		for(i in 0:target_sample){
			for(fi in 0:nj){
				x_down[i+1,fi+1] <-  SFS[fi+1] * choose(fi, i) * choose(nj - fi, target_sample - i) / choose(nj, target_sample)
			}
		}
		down_SFS <- rowSums(x_down)
		# down_SFS <- rowSums(ceiling(x_down))
		names(down_SFS) <- 0:target_sample
		if(!include_invariable) down_SFS <- down_SFS[-1]
		return(round(down_SFS))

	}
}



#' Get allele frequency from a filterVCF, VIF or VCFr object
#' 
#' If a filterVCF or VIF object is provided, invariant and fixed sites are automatically added to the frequency. Otherwise the user need to provide a Grange for both invariant and fixed sites. The function can use the ff library to reduce the allocated RAM used by the genotype matrix and the function extract_gt_ff() is used otherwise the more rapid extract.gt() from vcfR package is used.
#' 
#' The number of sample is taken from the vcfR object directly. the function returns a S3 list containing 3 GRanges element, snp frequecncy, invariant and fixed sites.
#' 
#' @param VCF a filterVCF object, a VIF object or a vcfR object
#' @param Invariant a Grange with invariant position on ly if VCF is of class "vcfR"
#' @param fixed a Grange with fixed position on ly if VCF is of class "vcfR"
#' @param hets_prop_allowed = NULL Threshold to remove sites with a certain proportion of heterozygous, intented for haploid.
#' @param use_ff = TRUE Should the genotype be ectracted using the ff library, this reduce the RAM memmory but is a bit slower
#' @return
#' \item{vcf_freq }{A Granges object for SNPs with frequencies as metadata} %% 
#' \item{invariant }{A Granges object for invariant with frequencies=0 as metadata} %%
#' \item{fixed }{A Granges object for fixed with frequencies=nb of chromosome sample in the population as metadata} %% 
#' ...
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords popstat
#' @examples
#' 
#' 
#' freq <- get_frequency(filterVCF_object)
#' sapply(freq, length)
#' 
#' @export
get_frequency <- function(VCF, Invariant = NULL , fixed = NULL, hets_prop_allowed = NULL, use_ff = TRUE){
	# invariant and fixed are file name
	if(class(VCF) == "filterVCF" | class(VCF) == "VIF"){
		Invariant <- invariant(VCF)
		fixed <- fixed(VCF)
		#add reference to the definition
		reference <- reference(VCF)
		VCF <- vcf(VCF)
	}

	if(isTRUE(use_ff)){
		gt <- extract_gt_ff(VCF)
		hom <- rowSums(sapply(1:ncol(gt), function(x) create_index_ff(gt[x], "==2"))) * 2
		het <- rowSums(sapply(1:ncol(gt), function(x) create_index_ff(gt[x], "==1")))
		nas <- rowSums(sapply(1:ncol(gt), function(x) create_index_ff(gt[x], "== -9")))
		freq <- hom + het
	}else{
		gt <- extract.gt(VCF, convertNA = TRUE, as.numeric = FALSE)
		hom 	<- rowSums(gt == "1/1" | gt == "1|1" , na.rm = TRUE) * 2
		het 	<- rowSums(gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0" , na.rm = TRUE)
		nas 	<- rowSums(is.na(gt), na.rm = TRUE)		
		freq <- hom + het
	}

	vcf_bed <- vcf2Grange(VCF)
	elementMetadata(vcf_bed)["freq"] <- freq
	elementMetadata(vcf_bed)["nsample"] <- (ncol(gt) - nas) * 2
	# add invariable
	# add fixed
	if(!is.null(hets_prop_allowed)){
		if(hets_prop_allowed < 1){
			hets_prop_allowed * ncol(gt)
		}
		index_het <- het > hets_prop_allowed
		vcf_bed <- vcf_bed[!index_het]
	}

	if(class(Invariant) != "GRanges"){
		Invariant <- bed2Grange(Invariant)
	}
	if(class(fixed) != "GRanges"){
		fixed <- bed2Grange(fixed)		
	}
	Invariant <- intersect_Granges(Invariant, fixed, invert = TRUE, merge_interval = FALSE)
	return(list(vcf_freq = vcf_bed, invariant = Invariant, fixed = fixed))
}



#' Compute pi , Wahterson theta, Tajima's D, Fu's Fs, Fay and Hu H and SFS in windows or per genomic ranges
#' 
#' This function call get_theta_pi_SFS() in windows or on each element of a Grange. Thus it can output stats for separate genes or any ranges defined in windows size. See get_theta_pi_SFS() for more info on the statistic computed
#' 
#' 
#' @param vcf_freq a vcf_freq object obtained qith get_frequency
#' @param reference A DNAStringSet with the reference genomes or the path to the fasta reference genome. If windows size is a GRanges set reference = NULL.
#' @param windows_size=20000 Size of the windows in bp. Alternatively A Granges with for wichi stats should be computed on each interval.
#' @param step = 20000 Size of the sliding in bp. Set step = windows_size for non-overlapping windows and ignore if windows_size is a GRange
#' @param downsample_size = 20 size of the downsampling of the SFS. Note that is only optional but the full SFS is always outputed.
#' @return A Granges object each range corresponding to either a window or a range from the inputed GRange (see above) with associated stats and information about the number of sites, total, segragating and fixed and the SFS

#' @author  ~~Benjamin Laenen~~
#' @references  
#' @examples
#' 
#' #Get the allele frequencies
#' freq <- get_frequency(filterVCF_object)
#' 
#' #import the gff and keep only CDS
#' gff_file <- "yourspecies.gff3"
#' gff <- import(gff_file, format = "gff")
#' gff <- gff[values(gff)[,"type"] == "CDS"] 
#' 
#' 
#' #Get the statistic for each CDS
#' CDS_stats <- get_stats_windows(freq, reference(filterVCF_object), windows_size = gff)
#' 
#' #Get the statistic in 50kb windows
#' windows_stats <- get_stats_windows(freq, reference(filterVCF_object), windows_size = 50000, step = 50000)
#' 
#' @export
get_stats_windows <- function(vcf_freq, reference, windows_size=20000, step = 20000, downsample_size = 20){
	# Compute pi , what theta, tajima D, Fu Fs, Fay and Hu, SFS in windows
	#
	if(class(windows_size) == "GRanges"){
		windows_genome <- windows_size
	}else{
		windows_genome <- windows_from_reference(reference, vcf_freq$vcf_freq, windows_size = windows_size)[[1]]
	}
    n <- max(elementMetadata(vcf_freq$vcf_freq)["nsample"][,1])

	index_vcf <- GenomicRanges::findOverlaps(vcf_freq$vcf_freq, windows_genome , select ="first",type="within", ignore.strand=TRUE)
	#vcf_freq_by_windows  <- split(vcf_freq$vcf_freq, index_vcf)
	vcf_freq_by_windows <- as.list(rep(NA, length(windows_genome)))
 	for(i in seq_along(windows_genome)) vcf_freq_by_windows[[i]] <- vcf_freq$vcf_freq[index_vcf == i & !is.na(index_vcf)]

	index_invariant <- GenomicRanges::findOverlaps(vcf_freq$invariant, windows_genome , select ="first", ignore.strand=TRUE)
	# invariant_by_windows <- split(vcf_freq$invariant, index_invariant)
	invariant_by_windows <- as.list(rep(NA, length(windows_genome)))
 	for(i in seq_along(windows_genome)) invariant_by_windows[[i]] <- vcf_freq$invariant[index_invariant == i & !is.na(index_invariant)]

	index_fixed <- GenomicRanges::findOverlaps(vcf_freq$fixed, windows_genome , select ="first", ignore.strand=TRUE)
	# fixed_by_windows     <- split(vcf_freq$fixed, index_fixed)
	fixed_by_windows <- as.list(rep(NA, length(windows_genome)))
 	for(i in seq_along(windows_genome)) fixed_by_windows[[i]] <- vcf_freq$fixed[index_fixed == i & !is.na(index_fixed)]

	stats_in_windows <- lapply(seq_along(windows_genome), function(i) get_theta_pi_SFS(vcf_freq_by_windows[[i]], invariant_by_windows[[i]], fixed_by_windows[[i]],  downsample_size = downsample_size, output_stats_only=TRUE, n = n))

	elementMetadata(windows_genome)["total_nb_sites"] <- unlist(sapply(stats_in_windows, "[", "total_nb_sites"))
	elementMetadata(windows_genome)["S"]              <- unlist(sapply(stats_in_windows, "[", "S"))
	elementMetadata(windows_genome)["REF_fixed"]      <- unlist(sapply(stats_in_windows, "[", "REF_fixed"))
	elementMetadata(windows_genome)["ALT_fixed"]      <- unlist(sapply(stats_in_windows, "[", "ALT_fixed"))
	elementMetadata(windows_genome)["nuc_div"]        <- unlist(sapply(stats_in_windows, "[", "nuc_div"))
	elementMetadata(windows_genome)["thetaW"]         <- unlist(sapply(stats_in_windows, "[", "thetaW"))
	elementMetadata(windows_genome)["tajimaD"]        <- unlist(sapply(stats_in_windows, "[", "tajimaD"))
	elementMetadata(windows_genome)["FayWu_H"]        <- unlist(sapply(stats_in_windows, "[", "FayWu_H"))
	elementMetadata(windows_genome)["Zeng_E"]         <- unlist(sapply(stats_in_windows, "[", "Zeng_E"))

	SFS_mat <- matrix(unlist(sapply(sapply(stats_in_windows, "[", "SFS"), "[", "SFS")), ncol = n+1, byrow = TRUE)
	SFS_name <- paste("SFS", 0 :n ,sep ="_")
	for(i in seq_along(SFS_name)){
		elementMetadata(windows_genome)[SFS_name[i]] <- SFS_mat[,i]
	}

	fSFS_mat <- matrix(unlist(sapply(sapply(stats_in_windows, "[", "SFS"), "[", "fSFS")), ncol = (n/2)+1, byrow = TRUE)
	fSFS_name <- paste("fSFS", 0 :(n/2) ,sep ="_")
	for(i in seq_along(fSFS_name)){
		elementMetadata(windows_genome)[fSFS_name[i]] <- fSFS_mat[,i]
	}

	down_SFS_mat <- matrix(unlist(sapply(sapply(stats_in_windows, "[", "SFS"), "[", "down_SFS")), ncol = downsample_size+1, byrow = TRUE)
	down_SFS_name <- paste("down_SFS", 0 :downsample_size ,sep ="_")
	for(i in seq_along(down_SFS_name)){
		elementMetadata(windows_genome)[down_SFS_name[i]] <- down_SFS_mat[,i]
	}

	down_fSFS_mat <- matrix(unlist(sapply(sapply(stats_in_windows, "[", "SFS"), "[", "down_fSFS")), ncol = downsample_size+1, byrow = TRUE)
	down_fSFS_name <- paste("down_fSFS", 0 :downsample_size ,sep ="_")
	for(i in seq_along(down_fSFS_name)){
		elementMetadata(windows_genome)[down_fSFS_name[i]] <- down_fSFS_mat[,i]
	}

	return(windows_genome)
}



#' Save output from Popstats into a BED file
#' 
#' More generally this function saves a Granges with metadata to a GRange and
#' is equivalent to Grange2bed(gr, output_name = "stat.bed", merge = FALSE,
#' keep_extra_col = TRUE)
#' 
#' @param gr a GRange or output from Popstat()
#' @param outputdir="." Output directory
#' @param output_name="stat.bed" name of the file
#' @return Save a BEDfile
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @examples
#' 
#' 
#' stats2bed(Popstat(filterVCFobject), windows_size = 20000)
#' 
#' @export
stats2bed <- function(gr, outputdir=".", output_name="stat.bed"){
	df <- data.frame(
		CHROM=seqnames(gr),
		START=start(gr)-1,
		END=end(gr),
		windows_nb=1:length(gr),
		strands=strand(gr))
    df <- cbind(df, elementMetadata(gr))

    write.table(df, file=paste(outputdir, output_name, sep ="/"), quote=F, sep="\t", row.names=F, col.names=TRUE)
}


#' Compute population genetic statistics and site frequency spectrum (SFS)
#' 
#' This function is used internally in Popstat() to compute, nucleotide
#' diversity, Whaterson's theta, Tajima's D, Fay and Wu H, and Zeng E. The SFS
#' is also extracted considering all sites that were sequenced including SNPs,
#' invariant and fixed and return the unfolded, folded and downsampled SFS. Note
#' that for unfolded SFS the user must first repolarize the vif object using
#' make_repolKey() and repolarize_VIF(). This function is used to calculate stats
#' on subset of the genome, such as gene, loci, windows and will always consider
#' the actual number of sites examined.
#' 
#' Nucleotide diversity Pi
#' 	Is calculated per sites using corrected Begun (2007) formula and divided by the number of sites actually sequenced. (Note : there was an error in the denominator in the original paper): 
#' 		pi <- 2 * (freq ALT * (nsample - freq ALT)) / (nsample * (nsample -1))
#' 	Begun, D. J., Holloway, A. K., Stevens, K., Hillier, L. W., Poh, Y. P., Hahn, M. W., ... & Pachter, L. (2007). Population genomics: whole-genome analysis of polymorphism and divergence in Drosophila simulans. PLoS biology, 5(11), e310.
#' 
#' Whaterson's theta
#' Is calculated as follow and divided by the number of sites actually sequenced:	
#' 		S / sum(1 / 1:(nsample - 1))
#' 	Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. Theoretical population biology, 7(2), 256-276.
#' 
#' Tajima's D
#' 	The numerator 'd' is calculated as d = sum(pi) - (S / sum(1 / 1:(nsample - 1))). 'd' is then passed too pegas::tajima.test() to get D. Note that for Tajima's D the total number of sites is irelevant as it does the difference between two estimator of theta this is unscaled.
#' 	Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
#' 
#' Fu's Fs 
#' 	Not yet implemented. (make a git merge request if you have some code)
#' 
#' ##NOT TESTED YET## used it with care	
#' Fay and Wu H
#' 	Calculated as follow:
#'		thetaH  <- sum(SFS[-c(1,length(SFS))] * (1:(nsample-1))^2) / choose(nsample,2)
#'		H <- sum(pi) - thetaH
#' Fay, J. C., & Wu, C. I. (2000). Hitchhiking under positive Darwinian selection. Genetics, 155(3), 1405-1413.
#' 
#' ##NOT TESTED YET## used it with care	
#'		thetaL  <- sum(SFS[-c(1,length(SFS))] * (1:(nsample-1))) / (nsample-1)^-1
#'		E <- thetaL - thetaW
#'	K. Zeng, Y.-X. Fu, S. Shi, and C.-I. Wu. Statistical tests for detecting positive selection by utilizing high-frequency variants. Genetics, 174:1431–1439, 2006.
#'
#'
#' @param vcf_freq a vcf_frequency object obtain using get_frequency() or GRanges for SNPs with frequency as value.
#' @param invariant GRange for invariant, not used if a vcf_freq object from get_frequency() is provided
#' @param fixed GRange for fixed, not used if a vcf_freq object from get_frequency() is provided
#' @param downsample_size = 20 Size of the downsampling
#' @param output_stats_only=TRUE return only the stats, if FALSE return a vcf_freq object and the stats
#' @param n=NULL Number of sample, if NULL it is taken from the vcf_freq. NOTE that this parameter influence all statistics, change it only if you know what you are doing!
#' @param SFS_only = FALSE report only the SFS
#' @return A list with statistics
#' \item{total_nb_sites }{Total number of sites sequenced, denominator of nuc_div and thetaW} 
#' \item{S }{Number of segragating sites} 
#' \item{REF_fixed }{Total number of invariant sites for the REF allele} 
#' \item{ALT_fixed }{Total number of fixed sites for the ALT allele} 
#' \item{SFS }{list containing  unfolded, folded and downsampled SFSs} 
#' \item{nuc_div }{Nucleotide diversity Pi relative to the total number of sites sequenced} 
#' \item{thetaW }{Whaterson's theta relative to the total number of sites sequenced} 
#' \item{tajimaD }{Tajima's D} 
#' \item{FayWu_H }{Fay and Wu H} 
#' \item{Zeng_E }{Zeng E} 
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~popstat
#' @examples
#' 
#' freq <- get_frequency(filterVCF_object)
#' get_theta_pi_SFS(freq, downsample_size = 10)
#' 
#' @export
get_theta_pi_SFS <- function(vcf_freq, invariant, fixed, downsample_size = 20, output_stats_only=TRUE, n=NULL, SFS_only = FALSE){
	if(all(names(vcf_freq) == c("vcf_freq", "invariant", "fixed")) & !is.null(names(vcf_freq))){
		invariant <- vcf_freq$invariant
		fixed <- vcf_freq$fixed
		vcf_freq <- vcf_freq$vcf_freq
	}
	if(length(vcf_freq) == 0){
		total_nb_sites=sum(width(invariant)) + sum(width(fixed)) + sum(width(vcf_freq))
		#return NA if there is no seq but 0 
		#if there is a seq but no variation 
		if(isTRUE(total_nb_sites == 0)){
			return(list(total_nb_sites=0,
						S=NA,
						REF_fixed=NA,
						ALT_fixed=NA,
						nuc_div=NA,
						thetaW=NA,
						tajimaD=NA,
						FayWu_H=NA,
						Zeng_E=NA,
						SFS=list(SFS=rep(NA, n+1),
								fSFS=rep(NA, (n/2)+1),
								down_SFS=rep(NA, downsample_size+1),
								down_fSFS=rep(NA, downsample_size+1))))
		}else{
			return(list(total_nb_sites=total_nb_sites,
						S=0,
						REF_fixed=sum(width(invariant)),
						ALT_fixed=sum(width(fixed)),
						nuc_div=0,
						thetaW=0,
						tajimaD=0,
						FayWu_H=0,
						Zeng_E=0,
						SFS=list(SFS=rep(0, n+1),
								fSFS=rep(0, (n/2)+1),
								down_SFS=rep(0, downsample_size+1),
								down_fSFS=rep(0, downsample_size+1))))
		}
	}else{
		stats <- list()
		total_nb_sites <- sum(width(invariant)) + sum(width(fixed)) + sum(width(vcf_freq))
		freq_alt <- as.numeric(elementMetadata(vcf_freq)["freq"][,1])
		nsample <- as.numeric(elementMetadata(vcf_freq)["nsample"][,1])
		if(is.null(n)) n = max(nsample)
		stats["total_nb_sites"] <- total_nb_sites
		stats["S"] <- sum(width(vcf_freq))
		stats["REF_fixed"] <- sum(width(invariant)) 
		stats["ALT_fixed"] <-  sum(width(fixed))

		## ==========================================================================
		#                          SFS
		# ==========================================================================

		possible_freq <- seq(0, n, by = 1)
		all_freq <- sort(unique(freq_alt))

		#Create a lookup table to assign frequency with NA to the right frequency
		lookup <- data.frame(find= all_freq, replace = NA)
		for(i in 1:length(all_freq)){
			lookup[i,"replace"] <- possible_freq[which.min(abs(possible_freq -  all_freq[i]))]
		}
		for(i in 1:nrow(lookup)){
			all_freq[all_freq == lookup[i, "find"]] <- lookup[i, "replace"]
		}

		#table of contengency counting each site in each frequency(), we are using tabulate instead of table to set the correct number of bins even if they are empty. Here the SFS must have  numbers of chromosome minus missing (they all have now the same number of chromosomes after the random adding of NA) and 1 extra bin for non-variant.
		SFS <- tabulate(freq_alt, nbins = n)

		#add no variation sites (a nb of sites)
		SFS <- c(sum(width(invariant)), SFS)

		#add the fixed alternates sites
		SFS[length(SFS)] <- sum(width(fixed))

		fSFS <- (SFS +rev(SFS))[1:ceiling(length(SFS)/2)]
		fSFS[ceiling(length(SFS)/2)] <- fSFS[ceiling(length(SFS)/2)] /2

		if(!is.null(downsample_size)){
			down_SFS <- downsample_SFS(SFS, downsample_size)
			down_fSFS <- downsample_SFS(fSFS, downsample_size)
		}else{
			down_SFS  <- 0
			down_fSFS <- 0
		}

		stats[["SFS"]] <- list(SFS=SFS, fSFS=fSFS, down_SFS=down_SFS, down_fSFS=down_fSFS)

		if(isTRUE(SFS_only)){
			return(stats)
		}
		# ==========================================================================
		#                     PI
		# ==========================================================================
		#check that carefully
		pi <- 2 * (freq_alt * (nsample - freq_alt)) / (nsample * (nsample -1))

		nuc_div <- sum(pi, na.rm=TRUE) / total_nb_sites

		stats["nuc_div"] <- nuc_div

		elementMetadata(vcf_freq)["pi"] <- pi

		# ==========================================================================
		#                     theta Wat
		# ==========================================================================

		#Site based Theta Watterson (1 / harm number n-1)
		harm <- function(n){sum(1 / 1:(n - 1))}

		elementMetadata(vcf_freq)["thetaW"] <- 1 / sapply(nsample, harm)
		#check if the sites is not invariant, fixed
		elementMetadata(vcf_freq)[nsample == freq_alt,"thetaW"] <- 0
		stats["thetaW"] <- sum(elementMetadata(vcf_freq)["thetaW"][,1]) / total_nb_sites

		# ==========================================================================
		#                     tajima D
		# ==========================================================================
		d = sum(pi) - stats[["S"]] / harm(n)
		S = sum(width(vcf_freq))
		stats["tajimaD"] <- tajimaD_test(d,n,S)$D

		# ==========================================================================
		#                     Fu Fs
		# ==========================================================================
		# library(strataG)
		# vcf_DNA <- vcfR2DNAbin(vcf, consensus = FALSE, extract.haps = TRUE, unphased_as_NA = FALSE, ref.seq = NULL,  verbose = FALSE)
		# fusFs(vcf_DNA)
		# vcf_gtype <- sequence2gtypes(vcf_DNA)
		# LDgenepop(vcf_gtype)

		# ==========================================================================
		#    Fay and Wu
		# ==========================================================================
		#Fay, J. C., & Wu, C. I. (2000). Hitchhiking under positive Darwinian selection. Genetics, 155(3), 1405-1413.
		##not sure yet
		thetaH  <- sum( SFS[-c(1,length(SFS))] * (1:(n-1))^2) / choose(n,2)
		H <- sum(pi) - thetaH
		stats["FayWu_H"] <- H
		# ==========================================================================
		#   Zeng et al.’s E
		# ==========================================================================
		# K. Zeng, Y.-X. Fu, S. Shi, and C.-I. Wu. Statistical tests for detecting positive selection by utilizing high-frequency variants. Genetics, 174:1431–1439, 2006.
		##experimental too
		thetaL  <- sum( SFS[-c(1,length(SFS))] * (1:(n-1))) / (n-1)^-1
		E <- thetaL - sum(elementMetadata(vcf_freq)["thetaW"][,1])
		stats["Zeng_E"] <- E

		# ==========================================================================
		if(output_stats_only){
			return(stats)
		}else{
			return(list(vcf_freq=vcf_freq, stats=stats))
		}

		}
}



#' Tajima's D test
#' 
#' Function accomodated from pegas. Used internally in get_theta_pi_SFS() and Popstat()
#' 
#'	Paradis, E. (2010). pegas: an R package for population genetics with an integrated–modular approach. Bioinformatics, 26(3), 419-420.
#' 
#' @param d numerator of Tajima's D : pi - thetaW
#' @param n number of chromosome sampled
#' @param S number of segragating sites
#' @return A list with three numeric values:
#' \item{D }{Tajima's D statistic.} 
#' \item{Pval.normal	 }{the p-value assuming that D follows a normal distribution with mean zero and variance one.} 
#' \item{Pval.beta }{the p-value assuming that D follows a beta distribution after rescaling on [0, 1] (Tajima, 1989).} 
#' 
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#'
#' @references  Tajima, F. (1989) Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123, 595–595.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' @export
tajimaD_test <- function(d, n, S){
	tmp <- 1:(n - 1)
	a1 <- sum(1/tmp)
	a2 <- sum(1/tmp^2)
	b1 <- (n + 1)/(3 * (n - 1))
	b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
	c1 <- b1 - 1/a1
	c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
	e1 <- c1/a1
	e2 <- c2/(a1^2 + a2)
	D <-  d /sqrt(e1 * S + e2 * S * (S - 1))
	Dmin <- (2/n - 1/a1)/sqrt(e2)
	Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
	tmp1 <- 1 + Dmin * Dmax
	tmp2 <- Dmax - Dmin
	a <- -tmp1 * Dmax/tmp2
	b <- tmp1 * Dmin/tmp2
	p <- pbeta((D - Dmin)/tmp2, b, a)
	if(!is.nan(p)){
		p <- if (p < 0.5) 2 * p else 2 * (1 - p)
	}
	list(D = D, Pval.normal = 2 * pnorm(-abs(D)), Pval.beta = p)
}

#' Create input for Sweepfinder2
#' 
#' This function create input for the program sweepfinder2 from a vcf_freq
#' object obtain using get_frequency()
#' 
#' 	# The input for sweepfinder looks like this
#'	# position x   n   folded
#'	# 460000   9   100 0
#'	# 460010   100 100 0
#'	# 460210   30  78  1
#'	# 463000   094 0
#'
#' 
#' @param vcf_freq a vcf_frequency object obtain using get_frequency()
#' @param output_file = "sweepfinder.input" <Name for the output file
#' @param invariant = TRUE include the invariant in the sweepfind2 input? The program can be run with or without.
#'
#' @return NULL
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~convert
#' @examples
#' 
#' freq <- get_frequency(filterVCF_object)
#' vcf_freq2sweepfinder(freq, output_file = "sweepfinder_all_sites.input") 
#' vcf_freq2sweepfinder(freq, output_file = "sweepfinder_snps.input",  invariant = FALSE) 
#' 
#' @export
vcf_freq2sweepfinder <- function(vcf_freq, output_file = "sweepfinder.input", invariant = TRUE){
	# The input for sweepfinder looks like this
	# position x   n   folded
	# 460000   9   100 0
	# 460010   100 100 0
	# 460210   30  78  1
	# 463000   094 0

	if(class(vcf_freq) == "filterVCF"){

	}
	max_freq <- max(elementMetadata(vcf_freq$vcf_freq)[,"nsample"])

	if(length(vcf_freq$fixed) != 0){
		vcf_freq$fixed <- expand_all_pos_Granges(vcf_freq$fixed)
		elementMetadata(vcf_freq$fixed)     <- data.frame(freq = rep(max_freq, length(vcf_freq$fixed)),     nsample = max_freq)
	}else{
		message("Fixed sites empty, check that is correct. Maybe they are included in the VCF?")	
	}
	
	if(isTRUE(invariant)){
		vcf_freq$invariant <- expand_all_pos_Granges(vcf_freq$invariant)
		elementMetadata(vcf_freq$invariant) <- data.frame(freq = rep(0,        length(vcf_freq$invariant)), nsample = max_freq)
		freq <- c(vcf_freq$vcf_freq, vcf_freq$invariant, vcf_freq$fixed)
	}else{
		freq <- c(vcf_freq$vcf_freq, vcf_freq$fixed)
	}
	freq <- sort(freq)

	sweepfinder_input <- data.frame(position = start(freq), 
									x = elementMetadata(freq)[,"freq"], 
									n = elementMetadata(freq)[,"nsample"], 
									folded = 1)
	write.table(sweepfinder_input, output_file, quote = FALSE, row.names = FALSE, sep = "\t")
}

#' Computed global population genetic statistics and window/range based
#' 
#' This function is the main utility to compute population genetic statsitics
#' using the filterVCF of VIF object class from this package. However, users can
#' also read in their own VCF file but must provide invariant sites as the main
#' purpose of this package is to always consider all sites that have bee
#' actually sequence when computing the statistic to avoid inflating value dues
#' to lack of sequencing in some genomic regions.
#' 
#' Global statistics will be calculated for the whole object. User can also
#' compute statisticsby windows based on the reference genome present in the
#' filterVCF or VIF object using the function get_stats_windows(). 
#' If a GRange is provided instead of window size, then
#' the statistics is calculated for each range (e.g. gene, exon, 4fold).
#' 
#' #' Nucleotide diversity Pi
#' 	Is calculated per sites using corrected Begun (2007) formula and divided by the number of sites actually sequenced. (Note : there was an error in the denominator in the original paper): 
#' 		pi <- 2 * (freq ALT * (nsample - freq ALT)) / (nsample * (nsample -1))
#' 	Begun, D. J., Holloway, A. K., Stevens, K., Hillier, L. W., Poh, Y. P., Hahn, M. W., ... & Pachter, L. (2007). Population genomics: whole-genome analysis of polymorphism and divergence in Drosophila simulans. PLoS biology, 5(11), e310.
#' 
#' Whaterson's theta
#' Is calculated as follow and divided by the number of sites actually sequenced:	
#' 		S / sum(1 / 1:(nsample - 1))
#' 	Watterson, G. A. (1975). On the number of segregating sites in genetical models without recombination. Theoretical population biology, 7(2), 256-276.
#' 
#' Tajima's D
#' 	The numerator 'd' is calculated as d = sum(pi) - (S / sum(1 / 1:(nsample - 1))). 'd' is then passed too pegas::tajima.test() to get D. Note that for Tajima's D the total number of sites is irelevant as it does the difference between two estimator of theta this is unscaled.
#' 	Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
#' 
#' Fu's Fs 
#' 	Not yet implemented. (make a git merge request if you have some code)
#' 
#' ##NOT TESTED YET## used it with care	
#' Fay and Wu H
#' 	Calculated as follow:
#'		thetaH  <- sum(SFS[-c(1,length(SFS))] * (1:(nsample-1))^2) / choose(nsample,2)
#'		H <- sum(pi) - thetaH
#' Fay, J. C., & Wu, C. I. (2000). Hitchhiking under positive Darwinian selection. Genetics, 155(3), 1405-1413.
#' 
#' ##NOT TESTED YET## used it with care	
#'		thetaL  <- sum(SFS[-c(1,length(SFS))] * (1:(nsample-1))) / (nsample-1)^-1
#'		E <- thetaL - thetaW
#'	K. Zeng, Y.-X. Fu, S. Shi, and C.-I. Wu. Statistical tests for detecting positive selection by utilizing high-frequency variants. Genetics, 174:1431–1439, 2006.
#' 
#' 
#' @param filterVCF A filterVCF or VIF object
#' @param vcf_file = NULL Unused if filterVCF is provided. path to a VCF file that can be read by read.vcfR()
#' @param invariant = NULL Unused if filterVCF is provided. A GRange with invariant sites.
#' @param fixed = NULL Unused if filterVCF is provided. A GRange with fixed sites.
#' @param reference = NULL Unused if filterVCF is provided. A DNAStringSet with the reference genome
#' @param downsample_size = NULL number of chromosome to downsample in the SFS. Note that the full SFS is always reported.
#' @param windows_size = NULL Size of the windows to compute statistics. If NULL only global stats is reported. If a GRange is provided, statistics will be calculated for each range (e.g. gene)
#' @param windows_step = NULL size of the sliding windows. If NULL non-overlapping windows is used
#' @param outputdir = NULL Name of the output directory if save2file = TRUE
#' @param output_name = NULL Root name of the output  if save2file = TRUE
#' @param sites = NULL A GRange to intersect with. For example, the user can provide non-synonymous sites and computes stats per windows or per genes.
#' @param save2file = FALSE Save output to a bed file with extra columns?
#' @param ... additional options from initialise_option()
#' @return A list with 3 elements
#' \item{vcf_freq }{A vcf_freq object}
#' \item{global_stats }{A list with the global statistics}
#' \item{stats_windows }{A Grange with stats in the metadata}
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~popstat
#' @examples
#' 
#' #Compute only global statistics
#' vif_object <- as.VIF(filterVCF_object)
#' Popstat(vif_object , windows_size = NULL, windows_step = NULL, outputdir = NULL, output_name = NULL, sites = NULL, save2file = FALSE)
#'  
#' #Compute only global statistics for synonymous
#' #Read in synonymous sites
#'  syn <- bed2Grange("4fold_syn.bed")
#' Popstat(vif_object , windows_size = NULL, windows_step = NULL, outputdir = NULL, output_name = NULL, sites = syn, save2file = FALSE)
#'  
#' #Compute statistics for synonymous per gene and save to a bed file
#' #Read in gff and selct only gene
#' gff <- import("my_species.gff3", format = "gff")
#' gff <- gff[values(gff)$type == "gene"]
#' 
#' Popstat(vif_object , windows_size = gff, windows_step = NULL, outputdir = "synonymous_stats", output_name = "my_species_syn_gene", sites = syn, save2file = TRUE)
#' 
#' #Compute statistics for all sites in 20kb windows and save to a bed file
#' Popstat(vif_object , windows_size = 20000, windows_step = NULL, outputdir = "all_sites_stats", output_name = "my_species_20kb", sites = NULL, save2file = TRUE)
#' 
#' # Use external VCF file from the user.
#' vcf_file = "my.vcf"
#' invariant <- bed2Grange("invariant_sites.bed")
#' #Note that it is not necesaary if they are still present in the VCF
#' fixed <- bed2Grange("fixed_sites.bed")
#' reference <- readDNAStringSet("my.genome.fasta")
#' 
#' Popstat(filterVCF = NULL ,  vcf_file = vcf_file, invariant = invariant, fixed = fixed, reference = reference, windows_size = 20000, windows_step = NULL, outputdir = "my_vcf_stats", output_name = "my_vcf_20kb", sites = NULL, save2file = TRUE)
#' 
#' @export
Popstat <- function(filterVCF = NULL, vcf_file = NULL, invariant = NULL, fixed = NULL, reference = NULL, downsample_size = NULL, windows_size = NULL, windows_step = NULL, outputdir = NULL, output_name = NULL, sites = NULL, save2file = FALSE,...){

	if(class(filterVCF) == "filterVCF" | class(filterVCF) == "VIF"){
		opt <- GetOpt(filterVCF)
		vcf <- vcf(filterVCF)
		invariant <- invariant(filterVCF)
		fixed <- fixed(filterVCF)
		#add reference to the definition
		reference <- reference(filterVCF)
		if(is.null(output_name)) output_name <- basename(ifelse(is.null(opt$vcf_file), "outStats", GetOpt(filterVCF)$vcf_file))
	}else{
		vcf <- read.vcfR(vcf_file, limit = 1e+08, verbose = TRUE, convertNA = FALSE)
		vcf@fix[,"ID"] <- NA
		if(class(invariant) != "GRanges" & file.exists(invariant)) invariant <- bed2Grange(invariant)
		if(class(fixed) != "GRanges" & file.exists(fixed)) fixed <- bed2Grange(fixed)
		if(is.null(output_name)) output_name <- basename(vcf_file)
	}

	if(!is.na(opt$chr[1])){
		vcf <- vcf[getCHROM(vcf) %in% opt$chr]
	}
		
	if(!is.na(opt$sample)){
		opt$sample_list <- parse_opt_sample(opt$sample)
		vcf <- vcf[,c("FORMAT", opt$sample_list)]
		warnings("Information in the INFO field will not be recalculated, use GATK -T SelectVariants on the output vcf to recreate them.")
	}

	if(!is.null(sites)){
		vcf <- intersect_Granges(vcf, sites)
		invariant <- intersect_Granges(invariant, sites)
		fixed <- intersect_Granges(fixed, sites) 
	}

	vcf_freq <- get_frequency(vcf, invariant, fixed)
	global_stats <- get_theta_pi_SFS(vcf_freq$vcf_freq, vcf_freq$invariant, vcf_freq$fixed, downsample_size, output_stats_only = FALSE, n = NULL)

	#updating vcf_freq with stats
	vcf_freq$vcf_freq <- global_stats$vcf_freq

	if(!is.null(windows_size)){
		stats_windows <- get_stats_windows(vcf_freq, reference, windows_size=windows_size, step = windows_step, downsample_size = downsample_size)
		if(isTRUE(save2file)){
			if(is.numeric(windows_size)) windows_size <- round(windows_size/1000,1) else windows_size <- "region"
			stats2bed(stats_windows, outputdir=outputdir, output_name=paste0(basename(output_name),"_", windows_size, "kb_", "stats.bed"))
		}
	}else{
		stats_windows <- NULL
	}
 	
	return(list(vcf_freq = vcf_freq, global_stats = global_stats, stats_windows = stats_windows))
}


