# ==========================================================================
#                               Converters
# ==========================================================================


# ' VIF to popGenome
# ' 
# ' WORK in progress
# '
# '
# ' @param VIF a VIF object
# ' @param vcffile 
# ' @return %% ~Describe the value returned %% If it is a LIST, use %%
# ' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
# ' 'comp2'} %% ...
# '
# ' @author  ~~Benjamin Laenen~~
# '
# ' @references  put references to the literature/web site here
# ' @keywords ~kwd1 ~kwd2
# ' @examples
# '
# # #' @export
# VIF2PopGenome <- function(VIF,vcffile = "vcf", ...){
# 	# Formal class 'GENOME' [package "PopGenome"] with 118 slots
# 	# GENOME.class@ BIG.BIAL             : list()
# 	# GENOME.class@ SLIDE.POS            : list()
# 	# GENOME.class@ big.data             : logi(0) 
# 	# GENOME.class@ gff.info             : logi(0) 
# 	# GENOME.class@ snp.data             : logi(0) 
# 	# GENOME.class@ basepath             : chr(0) 
# 	# GENOME.class@ project              : chr(0) 
# 	# GENOME.class@ populations          : list()
# 	# GENOME.class@ poppairs             : logi(0) 
# 	# GENOME.class@ outgroup             : logi(0) 
# 	# GENOME.class@ region.names         : chr(0) 
# 	# GENOME.class@ feature.names        : chr(0) 
# 	# GENOME.class@ genelength           : num(0) 
# 	# GENOME.class@ n.sites              : num(0) 
# 	# GENOME.class@ n.sites2             : num(0) 
# 	# GENOME.class@ n.biallelic.sites    : num(0) 
# 	# GENOME.class@ n.gaps               : num(0) 
# 	# GENOME.class@ n.unknowns           : num(0) 
# 	# GENOME.class@ n.valid.sites        : num(0) 
# 	# GENOME.class@ n.polyallelic.sites  : num(0) 
# 	# GENOME.class@ trans.transv.ratio   : num(0) 
# 	# GENOME.class@ keep.start.pos       : num(0) 
# 	# GENOME.class@ Coding.region        : num(0) 
# 	# GENOME.class@ UTR.region           : num(0) 
# 	# GENOME.class@ Intron.region        : num(0) 
# 	# GENOME.class@ Exon.region          : num(0) 
# 	# GENOME.class@ Gene.region          : num(0) 
# 	# GENOME.class@ Pop_Neutrality       : list()
# 	# GENOME.class@ Pop_FSTN             : list()
# 	# GENOME.class@ Pop_FSTH             : list()
# 	# GENOME.class@ Pop_Linkage          : list()
# 	# GENOME.class@ Pop_Slide            : list()
# 	# GENOME.class@ Pop_MK               : list()
# 	# GENOME.class@ Pop_Detail           : list()
# 	# GENOME.class@ Pop_Recomb           : list()
# 	# GENOME.class@ Pop_Sweeps           : list()
# 	# GENOME.class@ FSTNLISTE            : list()
# 	# GENOME.class@ nucleotide.F_ST      : num[0 , 0 ] 
# 	# GENOME.class@ nucleotide.F_ST2     : num[0 , 0 ] 
# 	# GENOME.class@ nuc.diversity.between: num[0 , 0 ] 
# 	# GENOME.class@ nuc.diversity.within : num[0 , 0 ] 
# 	# GENOME.class@ nuc.F_ST.pairwise    : num[0 , 0 ] 
# 	# GENOME.class@ nuc.F_ST.vs.all      : num[0 , 0 ] 
# 	# GENOME.class@ n.haplotypes         : num[0 , 0 ] 
# 	# GENOME.class@ hap.diversity.within : num[0 , 0 ] 
# 	# GENOME.class@ hap.diversity.between: num[0 , 0 ] 
# 	# GENOME.class@ Pi                   : num[0 , 0 ] 
# 	# GENOME.class@ PIA_nei              : num[0 , 0 ] 
# 	# GENOME.class@ haplotype.counts     : num[0 , 0 ] 
# 	# GENOME.class@ haplotype.F_ST       : num[0 , 0 ] 
# 	# GENOME.class@ hap.F_ST.pairwise    : num[0 , 0 ] 
# 	# GENOME.class@ Nei.G_ST.pairwise    : num[0 , 0 ] 
# 	# GENOME.class@ hap.F_ST.vs.all      : num[0 , 0 ] 
# 	# GENOME.class@ Nei.G_ST             : num[0 , 0 ] 
# 	# GENOME.class@ Hudson.G_ST          : num[0 , 0 ] 
# 	# GENOME.class@ Hudson.H_ST          : num[0 , 0 ] 
# 	# GENOME.class@ Hudson.K_ST          : num[0 , 0 ] 
# 	# GENOME.class@ Hudson.Snn           : num[0 , 0 ] 
# 	# GENOME.class@ Phi_ST               : num[0 , 0 ] 
# 	# GENOME.class@ hap.pair.F_ST        : num[0 , 0 ] 
# 	# GENOME.class@ MKT                  : num[0 , 0 ] 
# 	# GENOME.class@ Tajima.D             : num[0 , 0 ] 
# 	# GENOME.class@ SLIDE                : num[0 , 0 ] 
# 	# GENOME.class@ Fay.Wu.H             : num[0 , 0 ] 
# 	# GENOME.class@ Zeng.E               : num[0 , 0 ] 
# 	# GENOME.class@ theta_Tajima         : num[0 , 0 ] 
# 	# GENOME.class@ theta_Watterson      : num[0 , 0 ] 
# 	# GENOME.class@ theta_Fu.Li          : num[0 , 0 ] 
# 	# GENOME.class@ theta_Achaz.Watterson: num[0 , 0 ] 
# 	# GENOME.class@ theta_Achaz.Tajima   : num[0 , 0 ] 
# 	# GENOME.class@ theta_Fay.Wu         : num[0 , 0 ] 
# 	# GENOME.class@ theta_Zeng           : num[0 , 0 ] 
# 	# GENOME.class@ Fu.Li.F              : num[0 , 0 ] 
# 	# GENOME.class@ Fu.Li.D              : num[0 , 0 ] 
# 	# GENOME.class@ Yach                 : num[0 , 0 ] 
# 	# GENOME.class@ n.segregating.sites  : num[0 , 0 ] 
# 	# GENOME.class@ Rozas.R_2            : num[0 , 0 ] 
# 	# GENOME.class@ Fu.F_S               : num[0 , 0 ] 
# 	# GENOME.class@ Strobeck.S           : num[0 , 0 ] 
# 	# GENOME.class@ Kelly.Z_nS           : num[0 , 0 ] 
# 	# GENOME.class@ Rozas.ZZ             : num[0 , 0 ] 
# 	# GENOME.class@ Rozas.ZA             : num[0 , 0 ] 
# 	# GENOME.class@ Wall.B               : num[0 , 0 ] 
# 	# GENOME.class@ Wall.Q               : num[0 , 0 ] 
# 	# GENOME.class@ mult.Linkage         : num[0 , 0 ] 
# 	# GENOME.class@ RM                   : num[0 , 0 ] 
# 	# GENOME.class@ CL                   : num[0 , 0 ] 
# 	# GENOME.class@ CLmax                : num[0 , 0 ] 
# 	# GENOME.class@ CLR                  : num[0 , 0 ] 
# 	# GENOME.class@ MDSD                 : num[0 , 0 ] 
# 	# GENOME.class@ MDG1                 : num[0 , 0 ] 
# 	# GENOME.class@ MDG2                 : num[0 , 0 ] 
# 	# GENOME.class@ D                    : num[0 , 0 ] 
# 	# GENOME.class@ BD                   : num[0 , 0 ] 
# 	# GENOME.class@ BDF                  : num[0 , 0 ] 
# 	# GENOME.class@ BDF_bayes            : num[0 , 0 ] 
# 	# GENOME.class@ alpha_ABBA           : num[0 , 0 ] 
# 	# GENOME.class@ alpha_BABA           : num[0 , 0 ] 
# 	# GENOME.class@ beta_BBAA            : num[0 , 0 ] 
# 	# GENOME.class@Bd_clr
# 	# GENOME.class@Bd_dir
# 	# GENOME.class@P.Bd_clr
# 	# GENOME.class@f
# 	# GENOME.class@RNDmin
# 	# GENOME.class@D.z
# 	# GENOME.class@D.pval
# 	# GENOME.class@BDF.z
# 	# GENOME.class@BDF.pval
# 	# GENOME.class@BDF.SE
# 	# GENOME.class@D.SE
# 	# GENOME.class@jack.knife
# 	# GENOME.class@missing.freqs
# 	# GENOME.class@n.fixed.sites
# 	# GENOME.class@n.shared.sites
# 	# GENOME.class@n.monomorphic.sites
# 	# GENOME.class@genes
# 	# GENOME.class@region.data
# 	# GENOME.class@region.stats
	
# 	# GENOME.class@region.data@biallelic.matrix[[1]][,]
# 	# opening ff /private/var/folders/t8/qrz65n454f36dwxqcg3h3jrh0000gn/T/RtmpvFvPQ2/ff13f718a6e132.ff
# 	#           254 273 333 363 408 473 523 568 569 586 610 659 767 799 888 893 898 932 965 1018 1063 1419 1422 1787 1800 1941 1992 2131 3455 3464 3595 3655 3713 3796 3899 3901 3955 4110 4111 4232 4272 4391 4399 4431 4879 4973 4998 5041 5095 5200 5273 5276 5360 5424 5593 5723 5736 5798 5811 5835 5879 5911 5926
# 	# 105         0   0   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    1    0    0    0    0    0    0    0
# 	# 105.2       0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    1    0    0    0    0    1    0    0    0    0    0    0    1
# 	# 11          0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    1    0    0    0    0    0    1    1    0    0    0    0    0    0
# 	# 11.2        0   0   0   1   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    1    1    1    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    1    0    0    0    0    0    0    0    0    0
# 	# 1198        0   0   1   0   1   0   0   0   0   0   0   0   1   0   0   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0    0    0    0    0    0    0    1    1    0    0    0    0    0    0    0

# 	# bgzvcffile = paste0(vcffile, ".bgz")
# 	# gzvcffile = paste0(vcffile, ".gz")
# 	# bgzip(vcffile, overwrite = TRUE)
# 	# indexTabix(bgzvcffile, format ="vcf4")
# 	# GENOME.class <-  readVCF(bgzvcfile,numcols = 1, tid = "chr2", from = 1, topos = 4, include.unknown = TRUE)
# 	# GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE) 
# 	# GENOME.class <- diversity.stats(GENOME.class) 

# 	write.vcf(vcf(VIF), gzvcffile)
# 	gunzip(gzvcffile)
# 	outdir <- paste0("PopGenome_", vcffile)
# 	VCF_split_into_scaffolds(vcffile, outdir)
# 	GENOME.class<-readData(outdir, format="VCF", include.unknown=TRUE, ...)

# 	# GENOME.class <- new("GENOME")
# 	# lengthvcf <- 10
# 	# a <- greek_gene_pollen@vcf_filtered[1:lengthvcf,]
# 	# gt <- extract.gt(a)
# 	# pos <- getPOS(a)
# 	# sample_name <- colnames(gt)
# 	# sample_name <- paste0(sample_name, rep(c("", ".2"), length(sample_name)))
# 	# b <- strsplit(gt, "/|\\|")

# 	# index <- sapply(b, function(x) isTRUE(is.na(x)))
# 	# b[index] <- list(c(NaN, NaN))
# 	# biallelic.matrix <- matrix(t(as.numeric(unlist(b))), nrow = ncol(gt)*2, ncol = nrow(gt), dimnames = list(sample_name, pos))
# 	# if(GENOME.class@big.data) biallelic.matrix <- as.ff(biallelic.matrix)

# 	# GENOME.class@region.data@biallelic.matrix[[1]] <- biallelic.matrix

# 	#important to include invariant here!
# 	GENOME.class@n.sites <- sum(width(vcf(VIF))) + sum(width(fixed(VIF))) + nrow(vcf(VIF))

# 	# GENOME.class@region.data@populations
# 	# GENOME.class@region.data@populations2
# 	# GENOME.class@region.data@popmissing
# 	# #contains the row identifiers (biallelic.matrix) of the outgroup individuals
# 	# GENOME.class@region.data@outgroup
# 	# GENOME.class@region.data@outgroup2
# 	# GENOME.class@region.data@reading.frame
# 	# GENOME.class@region.data@rev.strand
# 	# GENOME.class@region.data@Coding.matrix
# 	# GENOME.class@region.data@Coding.matrix2
# 	# GENOME.class@region.data@UTR.matrix
# 	# GENOME.class@region.data@Intron.matrix
# 	# GENOME.class@region.data@Exon.matrix
# 	# GENOME.class@region.data@Gene.matrix
# 	# GENOME.class@region.data@CodingSNPS
# 	# GENOME.class@region.data@UTRSNPS
# 	# GENOME.class@region.data@IntronSNPS
# 	# GENOME.class@region.data@ExonSNPS
# 	# GENOME.class@region.data@GeneSNPS
# 	# GENOME.class@region.data@transitions
# 	# GENOME.class@region.data@biallelic.matrix
# 	# GENOME.class@region.data@biallelic.sites
# 	# GENOME.class@region.data@biallelic.sites2
# 	# GENOME.class@region.data@reference
# 	# GENOME.class@region.data@matrix_codonpos
# 	#include information from bed file vector of length=n.snps. TRUE:synonymous, FALSE:non-synonymous,NA:non-coding re- gion
# 	# GENOME.class@region.data@synonymous
# 	# GENOME.class@region.data@matrix_freq
# 	# GENOME.class@region.data@n.singletons
# 	# GENOME.class@region.data@trans.transv.ratio
# 	# GENOME.class@region.data@n.valid.sites
# 	# GENOME.class@region.data@n.biallelic.sites
# 	# GENOME.class@region.data@polyallelic.sites
# 	GENOME.class@region.data@n.nucleotides <- list(rep(GENOME.class@n.sites, ncol(VIF@vcf_filtered@gt)-1))
# 	# GENOME.class@region.data@biallelic.compositions
# 	# GENOME.class@region.data@biallelic.substitutions
# 	# GENOME.class@region.data@minor.alleles
# 	# GENOME.class@region.data@codons
# 	# GENOME.class@region.data@sites.with.gaps
# 	# GENOME.class@region.data@sites.with.unknowns
# 	# GENOME.class@region.data@included
# 	get.sum.data(GENOME.class) 
# 		#1=C 2=T 3=G 4=A
# 	return(GENOME.class)

# }

#' Convert a vcfR to plink PED/MAP format
#' 
#' vcfR to plink PED/MAP
#' 
#' WORK in progress, check the output.
#'
#' @param vcfR an vcfR object
#' @param output "plink"
#' @param population dataframe with samples and population assignment
#' @param gt optional genotype matrix extraxted with extract_gt_ff()
#' @return list
#' \item{ped }{character vector plink formated ped}
#' \item{map }{character vector plink formated map}
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @examples
#'
#' @export
vcfR2Plink <- function(vcfR, output = "plink", population = NULL, gt = NULL){
	samples <- colnames(vcfR@gt)[-1]
	if(is.null(population)){
		fam.id = "pop1"
	}else{
		population <- population[population[,1] %in% samples,]
		rownames(population) <- population[,1]
		fam.id <- population[samples, 2]
	}
	chromosome = getCHROM(vcfR)
	position = getPOS(vcfR)
	allele.1 = getREF(vcfR)
	allele.2 = getALT(vcfR)
	
	#extract genotypes and replace by allelic base
	if(is.null(gt)){
		gt <- as.data.frame(extract_gt_ff(vcfR))
	}
	hom_ref <- paste(allele.1, allele.1, sep = " ")
	het <- paste(allele.1, allele.2, sep = " ")
	hom_alt <- paste(allele.2, allele.2, sep = " ")
	missing <- paste(0, 0, sep = " ")
	for(i in seq_len(ncol(gt))){
		gt[gt[,i] == 0, i] <- hom_ref[gt[,i] == 0]
		gt[gt[,i] == 1, i] <- het    [gt[,i] == 1]
		gt[gt[,i] == 2, i] <- hom_alt[gt[,i] == 2]
		gt[gt[,i] ==-9, i] <- missing
	}
	#chromosome should match human chrom or numeric
	#thus notation are 0, 1,2..99 X Y MT Z W XY
	if(any(grepl("chr|scaffold", unique(chromosome)))){
		message(sprintf("Chromosome name should be numeric or X,Y,Z,MT in plink, converting %s to %s", unique(chromosome), gsub("chr|scaffold", "", unique(chromosome))))
		chromosome <- gsub("chr|scaffold", "", chromosome)
	}

	#create a ped file according to 
	#http://www.gwaspi.org/?page_id=145
	#ped <- data.frame(fam.id = fam.id, sample.id = samples, paternal.id = 0, maternal.id = 0, sex = -9, affection = 0, t(gt))
	# write.table(ped, paste0(output, ".ped"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
	message(sprintf("Writing plink ped/map to %s\n", paste(paste0(output, c(".ped", ".map"), collapse = " "))))
	collapsed_gt <- apply(t(gt),1, paste, collapse = " ")
	ped <- paste(fam.id = fam.id, sample.id = samples, paternal.id = rep(0, length(sample)), maternal.id = rep(0, length(sample)), sex = rep(-9, length(sample)), affection = rep(0, length(sample)), collapsed_gt, sep = " ")
	writeLines(ped, paste0(output, ".ped"))
	
	#create the map file
	map <- data.frame(chr = chromosome, marker.id = paste(chromosome, seq_along(chromosome), sep = ":"), gen.dist = 0, physical.dist = position)
	write.table(map, paste0(output, ".map"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

	#create a ped file according to 
	#http://www.gwaspi.org/?page_id=145

	return(list(ped = ped, map = map))
}