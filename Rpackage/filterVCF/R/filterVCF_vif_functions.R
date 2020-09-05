#' Create index of bootstraped vif object
#' 
#' This function does not produce a bootstrap of the original vif object as it
#' would be heavy on memmory to copy the entire vcf,invariable, fixed and divergence part of
#' the vif object. Instead it creates index of position and bootstrap those
#' indexes. The latter can be then used with get_bootstrap_VIF to produce the
#' bootstrapped dataset.
#' 
#' The bootstrap is a resample with replacement and be conducted using
#' different block.
#' 
#' The user can provide a GRange object with blocks that will be intersect
#' with the vif object, boostrapping with the same see the snp, invariable and
#' fixed sites ' If a region was already define the bootstrap_by can be "region"
#' and the function look in the vif object using bootGRange(vif) to get the
#' GRange object. This is useful if you already bootstrapped and then remove
#' some samples. ' Finally if bootstrap_by is a numeric, it define the windows
#' size in bp to create blocks. The functions looks up the reference in
#' reference(vif) and then create non-overlapping windows on which it perfomrs
#' the bootstrap
#' 
#' 
#' @param vif a vif object
#' @param nb_bootstrap Number of bootstrap to create index
#' @param bootstrap_by Either a Granges, "region" or a numeric indicative of windows size in "bp"
#' @param seed The seed to perform the resampling, if NULL the seed is taken randomly and store in bootSeed(vif)
#' @return A vif object with the slot @bootstrap filled with :
#' 	 \item{bootstrap matrix}{a bootstrap matrix with the index of each block. The matrix has a dimension of nb_block X nb_bootstrap. It can be retrieved using bootstrapMat(vif)}
#'   \item{number of bootstrap}{retrieved using Nboot(vif)}
#'   \item{type of bootstrap}{retrieved with bootBy(vif)}
#'   \item{indexes}{the indexes for the vcf, invariable and fixed. If divergenceare present they al get boostrap with the same seeds making it comparable among bootstrap.}
#' @author  ~~Benjamin Laenen~~
#' @keywords vif
#' @examples
#' 
#' vif_bootstraped <- bootstrap_VIF(vif, nb_bootstrap = 200, bootstrap_by = 20000)
#' head(bootstrapMat(vif_bootstraped))
#' 
#' @export
bootstrap_VIF <- function(vif, nb_bootstrap = 200, bootstrap_by = "region", seed = NULL){
	set.seed(seed)
	if(is.null(seed)){
		seed <- sample(1:2^18, nb_bootstrap)
	}
	bootSeed(vif) <- seed
	#create index to bootstrap not copying
	# We can boostrap by a set a ranges 
	# that we reshuffle with replacement.
	#  This can be :
	#  	1) a GRanges objects with custom ranges
	#  	2) "genes" the Granges used to intersect 
	#  	and present in GetOpt(vif)$intersect_Granges
	#  	assuming that each range represent a gene (or loci)
	#  	to bootstrap
	#  	3) a number for the windows size, then windows
	#  	are reshuffled 
	
	if(class(bootstrap_by) == "GRanges"){
		Granges2bootstrap <- bootstrap_by
	}
	if(class(bootstrap_by) == "character"){
		Granges2bootstrap <- bootGRange(vif)
	}
	if(class(bootstrap_by) == "integer"){
		Granges2bootstrap <- windows_from_reference(reference(vif), vcf = vcf(vif), windows_size = bootstrap_by)[[1]]
	}

	index_range <- 1:length(Granges2bootstrap)
	bootstrap_matrix <- sapply(1:nb_bootstrap, function(x) {
		set.seed(seed[x])
		sample(index_range, length(index_range), replace = TRUE)
	})
	
	#Create lists (vcf,inv,fixed) of index of position 
	#for each ranges to be able to shuffle the
	#elements and retrieve the vcf,inv and fixed bootstraped
	#dataset quickly.
	Hits_vcf_Grange <- GenomicRanges::findOverlaps(vcf2Grange(vcf(vif)), Granges2bootstrap, ignore.strand=TRUE, select="all")
	all_ranges_vcf <- as.list(rep(0, length(Granges2bootstrap)))
	all_ranges_vcf[unique(subjectHits(Hits_vcf_Grange))]  <- split(queryHits(Hits_vcf_Grange), subjectHits(Hits_vcf_Grange))

	Hits_inv_Grange <- GenomicRanges::findOverlaps(invariant(vif), Granges2bootstrap, ignore.strand=TRUE, select="all")
	all_ranges_inv <- as.list(rep(0, length(Granges2bootstrap)))
	all_ranges_inv[unique(subjectHits(Hits_inv_Grange))]  <- split(queryHits(Hits_inv_Grange), subjectHits(Hits_inv_Grange))

	Hits_fixed_Grange <- GenomicRanges::findOverlaps(fixed(vif), Granges2bootstrap, ignore.strand=TRUE, select="all")
	all_ranges_fixed <- as.list(rep(0, length(Granges2bootstrap)))
	all_ranges_fixed[unique(subjectHits(Hits_fixed_Grange))]  <- split(queryHits(Hits_fixed_Grange), subjectHits(Hits_fixed_Grange))

	Hits_div_Grange <- GenomicRanges::findOverlaps(divergence(vif), Granges2bootstrap, ignore.strand=TRUE, select="all")
	all_ranges_div <- as.list(rep(0, length(Granges2bootstrap)))
	all_ranges_div[unique(subjectHits(Hits_div_Grange))]  <- split(queryHits(Hits_div_Grange), subjectHits(Hits_div_Grange))

	bootstrapMat(vif) <- bootstrap_matrix
	bootGRange(vif) <- Granges2bootstrap
	vif@bootstrap$all_ranges_vcf <- all_ranges_vcf
	vif@bootstrap$all_ranges_inv <- all_ranges_inv
	vif@bootstrap$all_ranges_fixed <- all_ranges_fixed
	vif@bootstrap$all_ranges_div <- all_ranges_div
	bootBy(vif) <- bootstrap_by
	Nboot(vif) <- nb_bootstrap
	return(vif)
}


#' Extract SFS, folded SFS and divergence from a VIF object
#' 
#' The VIF object contains information on SNP, invariant and fixed sites by
#' default. Divergence can be added in the form of a GRange object with single
#' position and corresponding ingroup and outgroup(s) base.
#' 
#' The way the divergence is calculated depends wether the VIF object was
#' repolarized or not. If an outgroup is available and the VIF has been
#' repolarized using repolarize_VIF(), then divergence can be directly taken from the
#' SFS. The user can choose to count divergent sites as either all DERIVED sites
#' (including segregating) or only fixed DERIVED sites (instead of fixed
#' ALTERNATE to the reference when the VIF is not polarized). If an outgroup is
#' present but the VIF is not repolarized, the divergence can be calculated by
#' counting the number of sites that differ between the ingroup and outgroup
#' (NOTE that for it requires 2 outgroups and that the base is the same between
#' OUT1 and OUT2, it is hardcoded for the moment). It is still possible to
#' remove the SNPs from the divergence using this methods. The second methods
#' has the advantage of considering the whole genome alignment entirely while
#' the repolarized can only asses sites that were sequenced and called (SNP,
#' invariant and fixed). However, using repolarization offers a way to estimate
#' the probability that a sites is ancestral or derived. 
#' 
#' The user is responsible to check that the VIF has been reporlarized and no check is
#' performed. Note that creating the repolKey does not change the VIF directly
#' but just create the mean to do so, the VIF needs to be repolarized in a
#' second step.
#' 
#' @param vif a VIF object
#' @param downsample_size=20 report also a downsample SFS using downsample_sfs()
#' @param divergence = FALSE Should divergence be reported
#' @param use_repol = FALSE Is the VIF repolarized and should divergence be taken from the SFS?
#' @param remove_polym_from_div = FALSE remove polymorphism from divergence
#' @param freq = NULL provide a frequency object constructed using get_frequency(). This is meant to speed the retrieving of the SFS as the main speed bottleneck is getting the allele frequencies.
#' @return A list
#' \item{SFS }{unfolded SFS} %% 
#' \item{fSFS }{folded SFS} %% 
#' \item{down_SFS }{downsampled unfolded SFS} %% 
#' \item{down_fSFS }{downsampled folded SFS} %% 
#' \item{D }{number of divergent sites between ingroup and outgroup(s)} %% 
#' \item{lD }{total number of divergence} %% 
#' ...
#'
#' @author  ~~Benjamin Laenen~~
#' @keywords vif
#' @examples
#' 
#' getSFS_VIF(VIF_object) 
#' 
#' #include the divergence from a bed file in the VIF
#' divergence(VIF_object) <- bed2Grange("divergence.bed")
#' getSFS_VIF(VIF_object, divergence = TRUE) 
#' 
#' 
#' #make the repolarized key
#' repolKey(VIF_object)<- make_repol_key(VIF_object, methods = "ML")
#' repol_VIF_object <- repolarize_VIF(VIF_object)
#' getSFS_VIF(repol_VIF_object, divergence = TRUE, use_repol = TRUE) 
#' 
#' 
#' @export
getSFS_VIF <- function(vif, downsample_size=20, divergence = FALSE, use_repol = FALSE, remove_polym_from_div = FALSE, freq = NULL){
	if(nrow(vcf(vif)) == 0 ){
		if(length(invariant(vif)) == 0){
			return(NULL)
		}else{
			SFS <- list(SFS = c(sum(width(invariant(vif))), sum(width(fixed(vif)))))
		}
	}else{
		#try to spped the code by providing the frequency
		#this is a huge bottleneck to recompute it each time
		#getSFS_VIF is called
		if(is.null(freq)){
			SFS <- get_theta_pi_SFS(get_frequency(vif, use_ff = FALSE), downsample_size = downsample_size, SFS_only = TRUE)$SFS
		}else{
			SFS <- get_theta_pi_SFS(freq, invariant = invariant(vif) , fixed = fixed(vif) , downsample_size = downsample_size, SFS_only = TRUE)$SFS
		}
	} 
	
	if(divergence){
		#if we have repolarized the vcf 
		#we can use the fixed ancestral and derived counts
		#as divergence as they are relative to the outgroup.
		#
		#If we did not repolarize we can take the divergence
		#that have been intersected with the Grange to extract
		#the number of divergence and the number of sites 
		#that were alined minus the sites we removed for filtering
		
		if(use_repol){
			if(remove_polym_from_div){
				# if we remove the polymorphism from the
				# divergence, we only consider the last column
				# of the SFS and the total number do not include
				# segragating sites
				SFS$D <- SFS$SFS[length(SFS$SFS)]
				SFS$lD <- sum(SFS$SFS[c(1, length(SFS$SFS))])
			}else{
				# If we include polymorphism in the divergence
				# we counts segragating sites and fixed sites
				# as divergent. 
				SFS$D <- sum(SFS$SFS[-1])
				SFS$lD <- sum(SFS$SFS)
			}
		}else{
			#Take only the divergence that we could assess
			if(remove_polym_from_div){
				div <- intersect_Granges(divergence(vif), reduce(c(invariant(vif), fixed(vif))), merge_interval = FALSE)
			}else{
				div <- divergence(vif)
			}
			if(length(div) == 0){
				SFS$D <- 0
				SFS$lD <- 0				
			}else{
				indexREF_OUT1 <- elementMetadata(div)[,1] != elementMetadata(div)[,2]
				#require the out1 to equal out2
				#I dont considered using only one or three
				#outgroups here. But we can introduce an if
				#to use it or not
				indexOUT1_OUT2 <- elementMetadata(div)[,2] == elementMetadata(div)[,3]
				SFS$D <- sum(indexREF_OUT1 & indexOUT1_OUT2)
				SFS$lD <- sum(width(div[indexOUT1_OUT2]))
			}
		}
	}
	#gc()
	return(SFS)

}

#' Get one bootstrap replicates of VIF
#' 
#' Bootstrap is not performed directly when using bootstrap_VIF() to avoid
#' overflowing the memmory with repeated objects. Instead index to create
#' bootstrap repliactes are created in order to get boostraped VIF on the fly.
#' This function retrieved one boostrap using the index created in the first
#' step.
#' 
#' 
#' @inheritParams getSFS_VIF
#' @param index index of a single bootstrap to extract, should be lower than NBoot(vif)
#' @return a new VIF object that has been bootstrapped. Note that is not a
#' proper VIF and will not be compatible with all functions (e.g. intersect)
#' because it contains duplicated lines as it is a bootstrap.
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @keywords ~bootstrap
#' @examples
#' 
#' #create the index for 20kb windows bootstrap
#' vif <- bootstrap_VIF(vif, nb_bootstrap = 200, bootstrap_by = 20000)
#' #get the first replicate
#' boot_vif <- list()
#' boot_vif[[1]]  get_bootstrap_VIF(vif, index = 1)
#' #extract SFS from boostrap replicate
#' get_SFS_vif(boot_vif[[1]])
#' 
#' @export
get_bootstrap_VIF <- function(vif, index = 1){
	# intersect the vcf with each range
	# separately. It will result in duplicate sites
	# and it is expected! Check down the line 
	# if it is ok
	# Merge all small vcf into one
	# This use a bit toom uch energy
	# as the VCF ust be copied each time
	#message(sprintf("Getting bootstrap nb %s", index))
	bootstraped_vif <- new("VIF")

	# Instead I use which to get the index of which sites
	# intersect with each ranges
	# and then I just have to select with []
	# which sites to keep (or duplicate) in the vcf
	INDEX <- bootstrapMat(vif)[,index]
	boot_all_ranges_vcf        <- vif@bootstrap$all_ranges_vcf[INDEX]
	vcf(bootstraped_vif)       <- vcf(vif)[unlist(boot_all_ranges_vcf)]
	rm(boot_all_ranges_vcf)
	
	#bootstrapped the invariant and the fixed
	boot_all_ranges_inv        <- vif@bootstrap$all_ranges_inv[INDEX]
	invariant(bootstraped_vif) <- invariant(vif)[unlist(boot_all_ranges_inv)]
	
	boot_all_ranges_fixed      <- vif@bootstrap$all_ranges_fixed[INDEX]
	fixed(bootstraped_vif)     <- fixed(vif)[unlist(boot_all_ranges_fixed)]
	
	#bootstrap the divergence
	boot_all_ranges_div        <- vif@bootstrap$all_ranges_div[INDEX]
	#elementMetadata(vif@divergence) <- NULL
	divergence(bootstraped_vif)<- divergence(vif)[unlist(boot_all_ranges_div)]

	reference(bootstraped_vif) <- reference(vif)
	alternate(bootstraped_vif) <- alternate(vif)
	repolKey(bootstraped_vif)  <- repolKey(vif)
	bootstraped_vif@bootstrap  <- vif@bootstrap
	GetOpt(bootstraped_vif)@opt <- GetOpt(vif)

	#gc()
	return(bootstraped_vif)
}

#' Get one bootsrap replicate of frequency object
#' 
#' This function is meant for internal speed up in make_input_DFE and should only be used by experienced user
#' 
#' 
#' @inheritParams getSFS_VIF
#' @inheritParams get_bootstrap_VIF
#' @return a boostraped frequency object
#'
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  put references to the literature/web site here
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' #create the index for 20kb windows bootstrap
#' vif <- bootstrap_VIF(vif, nb_bootstrap = 200, bootstrap_by = 20000)
#' freq <- get_frequency(vif)
#' #get the first replicate
#' boot_freq <- list()
#' boot_freq[[1]] <- get_bootstrap_freq(freq, vif, index = 1 )
#' 
#' @export
get_bootstrap_freq <- function(freq, vif, bootstrap_matrix = bootstrapMat(vif), index = 1){
	# intersect the vcf with each range
	# separately. It will result in duplicate sites
	# and it is expected! Check down the line 
	# if it is ok
	# Merge all small vcf into one
	# This use a bit toom uch energy
	# as the VCF ust be copied each time
	#message(sprintf("Getting bootstrap freq nb %s", index))
	INDEX <- bootstrap_matrix[,index]
	boot_all_ranges_vcf <- vif@bootstrap$all_ranges_vcf[INDEX]
	res                 <- freq$vcf_freq[unlist(boot_all_ranges_vcf)]
	rm(boot_all_ranges_vcf)

	return(res)
}

#' Create a formated input site frequency spectrum for different programs 
#' 
#' Only DFE-alplha and polyDFE are supported for the moment
#' 
#' 
#' @param sfs_sel a list for selected sites obtained using getSFS_VIF()
#' @param sfs_neut a list for neutral sites obtained using getSFS_VIF()
#' @param name The names to print as comment in the input
#' @param Dn = 0 number of selected divergent sites (not mandatory if the SFSs are lists obtained using getSFS_VIF())
#' @param LDn = 0 total number of selected divergence sites (not mandatory if the SFSs are lists obtained using getSFS_VIF())
#' @param Ds = 0 number of neutral divergent sites (not mandatory if the SFSs are lists obtained using getSFS_VIF())
#' @param NDs = 0 total number of neutral divergence sites (not mandatory if the SFSs are lists obtained using getSFS_VIF())
#' @param folded = FALSE should the SFS be folded or not. Note that the folding is not perfomed here, just the folded SFS is retrieved from the list.
#' @param downsampled = FALSE should the SFS be downsampled or not. Note that the downsampling is not perfomed here, just the downsampled SFS is retrieved from the list.
#' @param output = c("DFE-alpha", "polyDFE", "dadi", "fastsimcoal2") Type of output [only DFE-alpha and polyDFE are currently supported]
#' @return A character vector with the formated SFS input

#' @author  ~~Benjamin Laenen~~
#' @references  put references to the literature/web site here
#' @keywords sfs
#' @examples
#' 
#'	#get synonymous non synonymous sites
#'	syn <- bed2Grange("4fold.bed")
#'	nonsyn <- bed2Grange("0fold.bed")
#'
#'	#intersect with the VIF
#'	vif_syn <- intersect_vif_Grange(vif, syn)
#'	vif_nonsyn <- intersect_vif_Grange(vif, nonsyn)
#'
#'	#repolarize, note that the repolKey must be created before on the vif
#'	#and will be inherited when we intersect with syn and nonsyn
#'	vif <- repolarize_vif(vif, methods = "ML")
#'	vif_syn <- repolarize_vif(vif_syn, methods = "ML")
#'	vif_nonsyn <- repolarize_vif(vif_nonsyn, methods = "ML")
#'	
#'	#get the SFSs
#'	sfs_syn    <- getSFS_VIF(vif_syn,  , divergence = TRUE, use_repol = TRUE, remove_polym_from_div = TRUE)
#'	sfs_nonsyn <- getSFS_VIF(vif_nonsyn, divergence = TRUE, use_repol = TRUE, remove_polym_from_div = TRUE)
#'	#get input for DFE-alpha
#'	VIF_out <- sfs2DFE(sfs_nonsyn, sfs_syn, name = "input_gene_set1", folded = FALSE, downsampled = FALSE, output = "DFE-alpha")
#'
#' 
#' @export
sfs2DFE <- function(sfs_sel, sfs_neut, name = "VIF", Dn = 0, LDn = 0, Ds = 0, NDs = 0, folded = FALSE, downsampled = FALSE, output = c("DFE-alpha", "polyDFE", "dadi", "fastsimcoal2")){
	#name <- paste(gsub(".sfs", "", names(sfs_sel)), gsub(".sfs", "", names(sfs_neut)), sep = "_")
	
	if(class(sfs_sel) == "list"){
		if(!folded & !downsampled) sel <- sfs_sel$SFS
		if(!folded & downsampled) sel <- sfs_sel$down_SFS
		if(folded & !downsampled) sel <- c(sfs_sel$fSFS, rep(0, length(sfs_sel$fSFS)-1))
		if(folded & downsampled) sel <- c(sfs_sel$down_fSFS, rep(0, length(sfs_sel$down_fSFS)-1))
		Dn <- sfs_sel$D
		LDn <- sfs_sel$lD
	}else{
		sel <- sfs_sel
	}

	if(class(sfs_neut) == "list"){
		if(!folded & !downsampled) neut <- sfs_neut$SFS
		if(!folded & downsampled) neut <- sfs_neut$down_SFS
		if(folded & !downsampled) neut <- c(sfs_neut$fSFS, rep(0, length(sfs_neut$fSFS)-1))
		if(folded & downsampled) neut <- c(sfs_neut$down_fSFS, rep(0, length(sfs_neut$down_fSFS)-1))
		Ds <- sfs_neut$D
		LDs <- sfs_neut$lD
	}else{
		neut <- sfs_neut
	}

	if(class(sfs_sel) == "list" | class(sfs_neut) == "list") folded = FALSE

	if(output == "DFE-alpha"){
		div <- paste(c(LDn,LDs), c(Dn,Ds),sep= "\t", collapse="\n")
		nb_chrom <- paste0(length(sel)-1 , "\n")
		
		if(folded){
			fsel <- sel + rev(sel)
			fneut <- neut + rev(neut)
			fsel[ceiling(length(fsel) /2)] <- fsel[ceiling(length(fsel) /2)] / 2
			fneut[ceiling(length(fneut) /2)] <- fneut[ceiling(length(fneut) /2)] /2
			fsel[(ceiling((length(fsel) /2)) + 1) : length(fsel)] <- 0
			fneut[(ceiling((length(fneut)) /2) + 1) : length(fneut)] <- 0
			sfs <- paste0(paste(fsel, collapse = " "), "\n", paste(fneut, collapse = " "), "\n")
		}else{
			sfs <- paste0(paste(sel, collapse = " "), "\n", paste(neut, collapse = " "), "\n")
		}
		return(paste0(name,"\n", div, "\n", nb_chrom, sfs))
	}
	if(output == "polyDFE"){
		l1 <- paste0("#", name, "\n")
		l2 <- paste0("1\t1\t", length(sel)-1 , "\n")
		neut <- paste(as.numeric(c(neut[- c( 1)], sum(neut), Ds, LDs)), collapse = "\t")
		sel  <- paste(as.numeric(c(sel [- c( 1)],  sum(sel),  Dn, LDn)), collapse = "\t")
		
		return(paste0(l1,l2, neut, "\n", sel, "\n"))
	}
}


#' Wrapper to create input for DFE-alpha or polyDFE
#' 
#' This function is intented to automatically create input for the program
#' DFE-alpha or polyDFE using information in the VIF object and GRange object
#' for selected and neutral sites. If divergence is present in the VIF object,
#' repolarization can be performed using either parsimony or ML. Note that for
#' ML the function use the software of P. Keightley 'est-sfs' which should be
#' present in your $PATH and be callable from a system() call.
#' 
#' @inheritParams getSFS_VIF
#' @param nonsyn A Granges with selected sites, A path to a bed file with selected sites or a VIF object containing only selected sites
#' @param syn A Granges with neutral sites, A path to a bed file with neutral sites or a VIF object containing only neutral sites
#' @param name Root name to give in the input
#' @param output = c("DFE-alpha", "polyDFE") type of output
#' @param folded = FALSE fold the SFS?
#' @param downsample_size = NULL integer indicating the size of the downsampling
#' @param nb_bootstrap = NULL Number of desired bootstrap. If NULL, no bootstrap is performed (can be usefull for quicj debugging)
#' @param bootstrap_by="region" Should the bootstrap be perform by region given in bootGRange(vif) or by windows. If windows is desired an integer representing windows size in bp is required. 
#' @param mc.cores = 4 number of core to use.
#' @param use_repol = TRUE Should a repolarizion using the outgroup be performed?
#' @param remove_polym_from_div = TRUE Remove the contribution of polymorphism from the divergence
#' @param stat currently unused
#' @param bootstrapped_vif_syn = NULL advanced option, provided a list with boostrapped vif for neutral sites to improve perfomance when calling the function multiple times.
#' @param bootstrapped_vif_nonsyn = NULL advanced option, see above
#' @param bootstrapped_freq_syn = NULL advanced option, see above
#' @param bootstrapped_freq_nonsyn = NULL advanced option, see above
#' @return a large list with intersect element and output
#' \item{SFS }{unfolded SFS} %% 
#' \item{fSFS }{folded SFS} %% 
#' \item{down_SFS }{downsampled unfolded SFS} %% 
#' \item{down_fSFS }{downsampled folded SFS} %% 
#' \item{D }{number of divergent sites between ingroup and outgroup(s)} %% 
#' \item{lD }{total number of divergence} %% 

#' \item{vif_syn }{vif object for neutral sites} %% 
#' \item{vif_nonsyn }{vif object for selected sites} %% 
#' \item{VIF_out }{The actual input for DFE-alpha or polyDFE as a character vector} %% 
#' \item{sfs_syn }{SFS list for neutral sites} %% 
#' \item{sfs_nonsyn }{SFS list for selected sites} %% 
#' \item{bootstrapped_sfs_syn }{list of bootstrapped neutral SFS} %% 
#' \item{bootstrapped_sfs_nonsyn }{list of bootstrapped selected SFS} %% 
#' \item{bootstrapped_vif_syn }{ list of bootstrapped neutral VIF object} %% 
#' \item{bootstrapped_vif_nonsyn }{ list of bootstrapped selected VIF object} %% 
#' \item{bootstrapped_freq_syn }{ list of bootstrapped neutral frequency object} %% 
#' \item{bootstrapped_freq_nonsyn }{ list of bootstrapped selected frequency object} %% 
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords dfe
#' @examples
#'
#' @export
make_input_DFE <- function(vif, nonsyn, syn, name, output = c("DFE-alpha", "polyDFE"), folded = FALSE, downsample_size = NULL, nb_bootstrap = NULL, bootstrap_by="region", mc.cores = 4, use_repol = TRUE, remove_polym_from_div = TRUE, stat = TRUE, bootstrapped_vif_syn = NULL, bootstrapped_vif_nonsyn = NULL, bootstrapped_freq_syn = NULL, bootstrapped_freq_nonsyn = NULL){
	#bootstrap vif and div

	#get syn non syn sites
	if(class(syn) == "VIF"){
		vif_syn <- syn
	}else{
		vif_syn <- intersect_vif_Grange(vif, syn)
	}
	if(class(nonsyn) == "VIF"){
		vif_nonsyn <- nonsyn
	}else{
		vif_nonsyn <- intersect_vif_Grange(vif, nonsyn)
	}

	if(use_repol){
		vif <- repolarize_vif(vif, methods = opt$repolarize[1])
		vif_syn <- repolarize_vif(vif_syn, methods = opt$repolarize[1])
		vif_nonsyn <- repolarize_vif(vif_nonsyn, methods = opt$repolarize[1])
	}	
	#get sfs
	sfs_syn    <- getSFS_VIF(vif_syn,    downsample_size = downsample_size, divergence = TRUE, use_repol = use_repol, remove_polym_from_div = remove_polym_from_div)
	sfs_nonsyn <- getSFS_VIF(vif_nonsyn, downsample_size = downsample_size, divergence = TRUE, use_repol = use_repol, remove_polym_from_div = remove_polym_from_div)

	VIF_out <- sfs2DFE(sfs_nonsyn, sfs_syn, name = paste0(name, "_obs"), folded = folded, downsampled = !is.null(downsample_size), output = output)
		
	if(!is.null(nb_bootstrap)){
		vif <- bootstrap_VIF(vif, nb_bootstrap = nb_bootstrap, bootstrap_by = bootstrap_by)
		vif_syn    <- bootstrap_VIF(vif_syn   , nb_bootstrap = Nboot(vif), bootstrap_by = bootGRange(vif), seed = bootSeed(vif))
		vif_nonsyn <- bootstrap_VIF(vif_nonsyn, nb_bootstrap = Nboot(vif), bootstrap_by = bootGRange(vif), seed = bootSeed(vif))
		
		# bootstrapped_vif_syn <- mclapply(seq_len(nb_bootstrap), function(x) get_bootstrap_VIF(vif_syn, index = x), mc.cores = mc.cores)
		# bootstrapped_vif_nonsyn <- mclapply(seq_len(nb_bootstrap), function(x) get_bootstrap_VIF(vif_nonsyn, index = x), mc.cores = mc.cores)
		# very slow find a better solution!
		# foun:-D!!!
		if(is.null(bootstrapped_vif_syn)){
			rerun_syn <- FALSE
			bootstrapped_vif_syn    <- mclapply_try(seq_len(nb_bootstrap), function(x) get_bootstrap_VIF(vif_syn,    index = x), use_BiocParallel = TRUE, no_cores = 1)
		}else{
			rerun_syn <- TRUE
		}
		if(is.null(bootstrapped_vif_nonsyn)){
			rerun_nonsyn <- FALSE
			bootstrapped_vif_nonsyn <- mclapply_try(seq_len(nb_bootstrap), function(x) get_bootstrap_VIF(vif_nonsyn, index = x), use_BiocParallel = TRUE, no_cores = 1)
		}else{
			rerun_nonsyn <- TRUE
		}

		if(is.null(bootstrapped_freq_syn)){
			freq_syn <- get_frequency(vif_syn, use_ff = FALSE)
			bootstrapped_freq_syn    <- mclapply_try(seq_len(nb_bootstrap), function(x) get_bootstrap_freq(freq_syn, vif_syn, bootstrapMat(vif_syn),   index = x), use_BiocParallel = TRUE, no_cores = 1)
		}
		
		if(is.null(bootstrapped_freq_nonsyn)){
			freq_nonsyn <- get_frequency(vif_nonsyn, use_ff = FALSE)
			bootstrapped_freq_nonsyn    <- mclapply_try(seq_len(nb_bootstrap), function(x) get_bootstrap_freq(freq_nonsyn, vif_nonsyn, bootstrapMat(vif_nonsyn),   index = x), use_BiocParallel = TRUE, no_cores = 1)
		}

		bootstrapped_sfs_syn    <- mclapply_try(seq_len(nb_bootstrap), function(x) {
			if(is.null(bootstrapped_vif_syn[[x]])){
				message(sprintf("bootstrap %s failed return NULL", x))
				return(NULL)
			}else{
				#message(sprintf("Get sfs bootstrap %s", x))
				getSFS_VIF(bootstrapped_vif_syn[[x]],    downsample_size = downsample_size, divergence = TRUE, use_repol = use_repol, remove_polym_from_div = remove_polym_from_div, freq = bootstrapped_freq_syn[[x]])
			}
		}, use_BiocParallel = TRUE, no_cores = 1)


		bootstrapped_sfs_nonsyn    <- mclapply_try(seq_len(nb_bootstrap), function(x) {
			if(is.null(bootstrapped_vif_nonsyn[[x]])){
				message(sprintf("bootstrap %s failed return NULL", x))
				return(NULL)
			}else{
				#message(sprintf("Get sfs bootstrap %s", x))
				getSFS_VIF(bootstrapped_vif_nonsyn[[x]],    downsample_size = downsample_size, divergence = TRUE, use_repol = use_repol, remove_polym_from_div = remove_polym_from_div, freq = bootstrapped_freq_nonsyn[[x]])
			}
		}, use_BiocParallel = TRUE, no_cores = 1)
	
		#if there is no SNPs in the bootstrap then
		#the function getSFS_VIF report NULL
		#As we cant use those bootstrap we remove
		#from both syn and non_syn
		index_null_boot <- sapply(bootstrapped_sfs_syn, is.null) | sapply(bootstrapped_sfs_nonsyn, is.null)
		names(bootstrapped_sfs_syn) <- names(bootstrapped_sfs_nonsyn) <- seq_len(nb_bootstrap)
		bootstrapped_sfs_syn <- bootstrapped_sfs_syn[!index_null_boot]
		bootstrapped_sfs_nonsyn <- bootstrapped_sfs_nonsyn[!index_null_boot]
		
		boot_VIF_out <- list()
		for(i in seq_len(sum(!index_null_boot))){
			boot_VIF_out[[i]] <- sfs2DFE(bootstrapped_sfs_nonsyn[[i]], bootstrapped_sfs_syn[[i]] , name = paste0(name,"_", i), folded = folded, downsampled = !is.null(downsample_size), output = output)
		}
		VIF_out <- c(VIF_out, boot_VIF_out)
	}
	return(list(vif = vif, vif_syn = vif_syn, vif_nonsyn = vif_nonsyn, VIF_out = VIF_out, sfs_syn = sfs_syn, sfs_nonsyn = sfs_nonsyn, bootstrapped_sfs_syn = bootstrapped_sfs_syn, bootstrapped_sfs_nonsyn = bootstrapped_sfs_nonsyn, bootstrapped_vif_syn = if(!rerun_syn) bootstrapped_vif_syn, bootstrapped_vif_nonsyn = if(!rerun_nonsyn) bootstrapped_vif_nonsyn, bootstrapped_freq_syn = if(!rerun_syn) bootstrapped_freq_syn, bootstrapped_freq_nonsyn = if(!rerun_nonsyn) bootstrapped_freq_nonsyn))
}

