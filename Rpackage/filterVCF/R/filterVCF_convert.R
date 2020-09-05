# ==========================================================================
#             Convert format functions
# ==========================================================================


#' Tranform a bed file to a GRanges object
#' 
#' Read Bed file format and transform to genomic ranges GRanges object. Only
#' the three first column of the BED file are necessary, all extra column are
#' stored in a as S4Vectors accessible via elementMetadata(GRange_object)
#' 
#' 
#' @param bed_file The path to a bed file
#' @return GRanges
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @family utilities
#' @examples
#' 
#' my_sites <- bed2Grange("path_to_my_bed_file")
#' 
#' @export
bed2Grange <- function(bed_file, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}

	if(isTRUE(opt$verbose)) message(sprintf("Converting bed file %s to Grange", bed_file), appendLF=TRUE)
	
	bed <- tryCatch(read.table(bed_file, sep = determine_file_separator(bed_file), ...), error=function(e) NULL)

	if(is.null(bed)){
		if(isTRUE(opt$verbose)) message(sprintf("Empty file %s returning empty Grange", bed_file), appendLF=TRUE)
		bed_GRange <- GRanges()
	}else{
		colnames(bed)[1:3] <- c("CHROM", "START", "END")
		bed_GRange <- makeGRangesFromDataFrame(bed, starts.in.df.are.0based=TRUE, keep.extra.columns = TRUE)
		# Saving the rest of the bed file as a data.frame 
		# in metadata
	}
	return(bed_GRange)
}



#' Tranform a vcfR object to a GRanges object
#' 
#' Take a vcfR object and transform to genomic ranges GRanges object. ' Some
#' element can be store in the metadata of the GRanges as S4Vectors. If metadata
#' = "present" only the the pressence or missing of genotype is recorded. If
#' metadata = "genotype" is stored in a dataframe in the metadata.
#' 
#' 
#' @inheritParams create_filter_repeats_in_windows
#' @param metadata Should extra information be stored in the GRanges object. Default is NULL, only position are recorded. "present" to get the presence or missing of genotype. "genotype" to get the genotype matrix
#' @return GRanges
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @family utilities VCFutilities
#' @examples
#' 
#' vcfGrange <- vcf2Grange(vcf = vcfR)
#' vcfGrange_present <- vcf2Grange(vcf = vcfR, metadata = "present")
#' vcfGrange_genotype <- vcf2Grange(vcf = vcfR, metadata = "genotype")
#' 
#' elementMetadata(vcfGrange)
#' elementMetadata(vcfGrange_present)
#' elementMetadata(vcfGrange_genotype)
#' 
#' @export
vcf2Grange <- function(vcf, metadata = NULL, ...){
	# if(opt$verbose) message(sprintf("Converting vcf to Grange"), appendLF=TRUE)
	vcf_df_bed <- data.frame(CHROM = getCHROM(vcf), START = getPOS(vcf),  END = getPOS(vcf))
	if(nrow(vcf_df_bed) == 0){
		vcf_bed_GRange <- GRanges()
	}else{
		vcf_bed_GRange <- makeGRangesFromDataFrame(vcf_df_bed, starts.in.df.are.0based=FALSE)
	}

	if(!is.null(metadata)){
		if(metadata == "present"){
			meta <- !is.na(extract.gt(vcf))
			meta[meta] <- 1
			meta[!meta] <- 0
			elementMetadata(vcf_bed_GRange) <- meta
		}
		if(metadata == "genotype"){
			elementMetadata(vcf_bed_GRange) <- extract.gt(vcf)
		}
	}
	return(vcf_bed_GRange)
}

#' Extract genotype using memmory efficient ffdf data.frame
#' 
#' The genotype matrix exported by the vcfR package can be quite large in RAM
#' memmory and this function aimed to circumvent this issue by using the ff
#' package to use memory efficiently. The genotypes are extracted from the VCF
#' field GT wheter they are phased or not but phase information is not retained.
#' The genotype are coded as 
#' -9 = missing
#' 1 = heterozygous
#' 0 = homozygous for REF allele
#' 2 = homozygous for ALT allele
#' 
#' ~~ If necessary, more details than the description above ~~
#' 
#' @inheritParams create_filter_repeats_in_windows
#' @return a ffdf with genotype
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @references http://ff.r-forge.r-project.org/
#' @keywords ~utilities
#' @family utilities VCFutilities
#' @examples
#' 
#' 
#' gt <- extract_gt_ff(vcfR)
#' gt
#' 
#' @export
extract_gt_ff <- function(vcf, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}
	
	
	if(opt$verbose) message(sprintf("Extracting genotyping as ffdf"), appendLF=TRUE)
	gt <- vcf@gt[,-1]
	names_gt <- colnames(gt)

	gt[grepl("^0/1|^0\\|1|^1/0|^1\\|0", gt, perl = TRUE)] <- 1
	gt[grepl("^1/1|^1\\|1", gt, perl = TRUE)] <- 2
	gt[grepl("^0/0|^0\\|0", gt, perl = TRUE)] <- 0
	gt[grepl("^./.|.\\|.", gt, perl = TRUE)] <- -9
	gt <- as.ffdf(as.ff(as.byte(gt), dim = dim(gt)))
	colnames(gt) <- names_gt
	void <- invisible(close.ffdf(gt))
	gc(reset = TRUE)
	return(gt)
}



#' Save logical vector to BED file
#' 
#' This function write a BED file from a filter create with a filter function
#' that output a logical the same size as vcf. The GRanges from the vcf must be
#' provided as the filter does not have information on the chr and position but
#' match the one of the vcf. You can use the function vcf2Grange(vcfR) to
#' produce this file.
#' 
#' 
#' @param filter the logical vector created by a filter function
#' @param filter_name the name of the filter
#' @param vcf_bed_GRange A GRange object with the position of the VCF corresponding to the filter
#' @param outputdir = "." the output directory to save the file to
#' @param output_name the name of the BED file to save
#' @param ... <what param does>
#' @return Write to file
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @family utilities filters
#' @examples
#' 
#' vcf_bed_GRange <- vcf2Grange(vcfR)
#' gt <- extract_gt_ff(vcfR)
#' filter_for_missing <- create_filter_missing <- function(gt, allowed_missing_threshold = 0.2)
#' filter2bed(filter_for_missing, "missing_0.2", vcf_bed_GRange, output_name = "filter_missing_0.2.bed")
#' 
#' 
#' @export
filter2bed <- function(filter, filter_name, vcf_bed_GRange, outputdir = ".", output_name, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}
	
	if(!is.null(filter)){
		if(opt$verbose) message(sprintf("Save filter %s to %s", filter_name, paste0(outputdir, "/", output_name)), appendLF=TRUE)
		filter_Grange <- reduce(vcf_bed_GRange[filter])
		export(filter_Grange, paste0(outputdir, "/", output_name), "BED", ignore.strand = TRUE)
	}
}



#' Save GRanges to BED file
#' 
#' This function write a BED file from a GRanges object. It is a wrapper that
#' uses internally rtracklayer export to export to BED. Note that the intervals
#' are merged before saving to reduce file size. If you one to keep all interval
#' you should use the export from rtracklayer directly.
#' 
#' 
#' @param Grange_filter the logical vector created by a filter function
#' @param name the name of the filter
#' @param vcf_bed_GRange A GRange object with the position of the VCF corresponding to the filter
#' @inheritParams filter2bed
#' @return Write to file
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~utilities
#' @family utilities
#' @examples
#' 
#' # Exporting the position of the vcf in a BED format
#' vcf_bed_GRange <- vcf2Grange(vcfR)
#' Grange2bed(vcf_bed_GRange, "vcf", vcf_bed_GRange, output_name = "vcf.bed")
#' 
#' 
#' @export
Grange2bed <- function(Grange_filter, name= "", outputdir = ".", output_name, merge = TRUE, keep_extra_col = FALSE, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}
	output_file <- paste0(outputdir, "/", output_name)
	if(opt$verbose) message(sprintf("Save GRange : %s to \n%s", name, output_file), appendLF=TRUE)
	if(!is.null(Grange_filter)){
		if(isTRUE(merge)){
			rtracklayer::export(reduce(Grange_filter), output_file, "BED", ignore.strand = TRUE)
		}else{
			if(isTRUE(keep_extra_col)){
				df <- data.frame(CHROM=seqnames(Grange_filter),
								  START=start(Grange_filter)-1,
								  END=end(Grange_filter),
								  STRANDS=strand(Grange_filter))
				df <- cbind(df, values(Grange_filter))
				colnames(df)[1] <- paste0("#", colnames(df)[1])
				if(file.exists(output_file)) file.remove(output_file)
				writeLines(paste(colnames(df), collapse = "\t"), output_file)
				write.table(df, file=output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append = TRUE)
			}else{
				rtracklayer::export(Grange_filter, output_file, "BED", ignore.strand = TRUE)
			}
		}
	}
}



#' Transform GRange to filter
#' 
#' This function takes a GRanges object and the GRanges from the vcf to create
#' a logical vector matching the position of the vcf.
#' 
#' 
#' @param Grange a Genomic range object
#' @inheritParams filter2bed
#' @return a logical vector of the size of nrow(vcf)
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~utilities
#' @family utilities
#' @examples
#' 
#' vcf_bed_GRange <- vcf2Grange(vcfR)
#' filter <- Grange2filter(Grange, vcf_bed_GRange)
#' vcf[!filter]
#' 
#' @export
Grange2filter <- function(Grange, vcf_bed_GRange, ...){
	filter <- !is.na(GenomicRanges::findOverlaps(vcf_bed_GRange, Grange , select ="first", ignore.strand=TRUE))
	return(filter)
}

#' Create an index of a ff vector
#' 
#' ff vector cannot be directly compared using which or an expression. This
#' function takes an expression in the form of e.g. "==1| == -9" and apply it to
#' the ff vector.  This is used internally in filter functions and should not be
#' used directly unless the user have a good understanding of ff object.
#' 
#' 
#' @param ff a ff vector
#' @param expr a string with an expression to evaluate
#' @return logical vector TRUE for element for which the expression is TRUE
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' gt <- extract_gt_ff(vcfR)
#' #Create a filter for fix heterozygous in the sample
#' filter <- rowSums(sapply(1:ncol(gt), function(x) create_index_ff(gt[x], "==1| == -9"))) == ncol(gt)
#' 
#' 
#' @export
create_index_ff <- function(ff, expr){
	expr <- unlist(strsplit(expr, "\\|"))
	expr <- parse(text=paste(paste0("ff[]" , expr),collapse = "|"))
	index <- ffwhich(ff, eval(expr))[]
	index2 <- 1:nrow(ff)
	index2[index] <- 0
	index2[index2 != 0] <- 1
	return(abs(index2 -1))
}
