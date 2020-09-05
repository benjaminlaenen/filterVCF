# ==========================================================================
#                   Utilities
# ==========================================================================

#' Split a Genomic Ranges by another one
#' 
#' This function takes a GRanges object and find the overlap for each range in
#' a second GRanges creating a list with length equal to the number of ranges
#' in the second arguments. For the moment the first GRanges must list all sites
#' as the function is intended to intersect sites in ranges. The function 
#' expand_all_pos_Granges can be used to produce a GRanges with all sites listed
#' as single and is call internally if the unique(width(Grange1)) is different from
#' 1
#' 
#' For example, if the first GRanges are type of sites (e.g. synonymous,
#' stop-codons, ...) and the second GRanges has gene interval, the function will
#' bin the sites per genes and return a list with length equal to the number of
#' genes. If there is no overlap then an empty GRanges is return so to keep the
#' number of ranges in input and ouput equal. Note that there is no check for 
#' duplicates.
#' 
#' @param Grange1 A GRanges object that will be binned according to Grange2
#' @param Grange2 A GRanges object to intersect range in Grange1
#' @return A list of GRanges with length equal to length(Grange2) and containing
#' intersected element of Grange1
#'
#' @author  ~~Benjamin Laenen~~
#'  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' # Define the sites to bin
#' Grange1 <- GRanges(seqnames="chr1", ranges= 1:100)
#' #Define the interval for binning
#' Grange2 <- GRanges(seqnames="chr1", ranges= c(1:10, 60:70, 100:110, 120:140)) 
#' 
#' splitted_Grange1 <- split_Grange(Grange1, Grange2)
#' length(Grange2) == length(splitted_Grange1)
#' 
#' @export
split_Grange <- function(Grange1, Grange2){
	if(any(unique(width(Grange1)) != 1)){
		Grange1 <- expand_all_pos_Granges(Grange1)
	}
	split_factor <- findOverlaps(Grange1, Grange2, select ="first", ignore.strand=TRUE)
	Grange1_splitted <- split(Grange1, split_factor)
	
	#add the intervals of the Granges
	#that have no sites in the vcf.
	#The length of the results should be the length 
	#of the Grange
	missing_Grange <- seq_along(Grange2)[!seq_along(Grange2) %in% unique(split_factor)]
	empty_Grange_list <- lapply(seq_along(missing_Grange), function(x) GRanges())
	names(empty_Grange_list) <- missing_Grange

	all_splitted_Grange <- c(as.list(Grange1_splitted), empty_Grange_list)
	all_splitted_Grange <- all_splitted_Grange[order(as.numeric(names(all_splitted_Grange)))]
	return(all_splitted_Grange)
}


#' Remove overlaping ranges
#' 
#' Function to test and remove ranges that are 
#' overlapping by keeping the longest ranges and 
#' keeping the metadata.
#' 
#' This function is intended to remove overlapping gene
#' and keep the longest (primary) transcript. If a GFF is
#' read into a GRanges, it can remove all the exon/CDS parts but
#' keeping the gene definition.
#' 
#' @param Grange A GRanges object
#' @param ... argument ignore.strand can be passed to findOverlaps to ignore strandness inoverlap
#' @return A GRanges object with the smallest overlapping interval(s) removed
#'
#' @author  ~~Benjamin Laenen~~
#'  objects to See Also as \code{\link{help}}, 
#' @references  put references to the literature/web site here
#' @keywords ~utilities
#' @examples
#' 
#' Grange <- GRanges(seqnames="chr1", ranges= c(1:10, 1:4, 2:5, 1:2, 100:1000, 120:140)) 
#' Grange_simplified <- remove_overlapping_range(Grange)
#' 
#' Grange <- GRanges(seqnames="chr1", ranges= c(1:10, 1:4, 2:5, 1:2, 100:1000, 120:140), strand = c("+", "+", "-", "-", "+", "+")) 
#' Grange_simplified <- remove_overlapping_range(Grange)
#' # Ignoring strand
#' Grange_simplified <- remove_overlapping_range(Grange, ignore.strand = TRUE)
#' 
#' @export
remove_overlapping_range <- function(Grange, ...) {
	if(!isDisjoint(Grange)){
		message("Some ranges are overlapping, only the longest overlapping range will be kept!")
		Overlap=findOverlaps(Grange, ...)
		index = which(from(Overlap) != to(Overlap))
		# For each overllaping ranges check
		# which one is thelongest and report the other
		# one to remove from the list
		smaller_overlapping_ranges2remove <- unique(sapply(index,  function(i) {
			fromto <- which.max(c(width(Grange[from(Overlap)[i]]), width(Grange[to(Overlap)[i]])))
			if(fromto == 1){
				return(to(Overlap)[i])
			}
			if(fromto == 2){
				return(from(Overlap)[i])
			}
		}))
		if(ncol(elementMetadata(Grange)) != 0) Grange_name <- elementMetadata(Grange[smaller_overlapping_ranges2remove])[,1] else Grange_name <- NULL
		message(sprintf("The following ranges are overlapping and will be removed:\n%s", paste(paste(as.character(seqnames(Grange[smaller_overlapping_ranges2remove])), start(Grange[smaller_overlapping_ranges2remove]), end(Grange[smaller_overlapping_ranges2remove]), Grange_name), collapse = "\n")))

		Grange <- Grange[-smaller_overlapping_ranges2remove]
	
	}
	#recursively remove ranges
	#if did not work
	#could end up in an infinite loop!
	if(!isDisjoint(Grange)){
		Grange <- remove_overlapping_range(Grange)
	}
	return(Grange)
}


#' Escape space and parenthesis in a path
#' 
#' output correct path when there is space or
#' parenthesis in the original path as in 
#' default when using dropbox pro
#' 
#' 
#' @param path A character with the path
#' @return A character with escaped path
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' 
#' 
#' @export
escape_dropbox_pro_path <- function(path){
	path <- gsub(" ", "\\ ", path, fixed=T)
	path <- gsub("(", "\\(", path, fixed=T)
	path <- gsub(")", "\\)", path, fixed=T)
	return(path)
}



#' @export
dry_run <- function(...){
	#testing if files exist
	file.exists(opt$vcf_file)
	#...to be continued
}

#' @export
textbox <- function(text, wide = 80){
	text_size <- nchar(text)
	if(text_size > wide) text_size <- wide
	hash_line <- paste(rep("#", text_size+4), collapse = "")
	text2 <- paste(strwrap(text, width = wide, initial = "", prefix=""), collapse ="\n")
	
	text2 <- strsplit(text2, "\n")[[1]]
	if(length(text2) > 1){
		#space2add <- text_size - sapply(text2, nchar)
		space2add <- wide - sapply(text2, nchar)
		space2add[space2add < 0] <- 0
		text2 <- paste(sapply(seq_along(text2), function(i) paste("#", text2[i], paste(rep(" ", space2add[i]), collapse = ""), "#")), collapse = "\n")
		hash_line <- paste0(hash_line, "#")
	}else{
		text2 <- paste0("# ", text, " #")
	}
	
	return(paste("", hash_line, text2, hash_line, "", "", sep = "\n"))
}


#' Get the call from Rscript
#' 
#' This function is not intended for the end-user
#' in the console. Print the call when using Rscript
#' 
#' ~~ If necessary, more details than the description above ~~
#' 
#' @return A character of length 1
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' 
#' @export
get_call <- function(){
	args <- commandArgs(trailingOnly=F)
	script_name <- gsub("--file=", "", grep("--file", args, value = TRUE))
    args <- R.utils::commandArgs(trailingOnly=T, asValues = TRUE)

	for(i in seq_along(args)){
		args[[i]] <- ifelse(nchar(names(args)[i]) == 1 , 
								ifelse(is.logical(args[[i]]), 
									paste0("-", names(args)[i]) , 
									paste0("-", names(args)[i], " '", args[[i]], "'")), 
								ifelse(is.logical(args[[i]]), 
									paste0("--", names(args)[i]), 
									paste0("--", names(args)[i], " '", args[[i]], "'"))
								)
	}

	call <- paste(c("Rscript", script_name, unlist(args)), collapse = " ")

	return(call)
}




#' Make windows from a reference
#' 
#' This function creates non-overlapping or slidding windows from a reference. 
#' 
#' 
#' @param reference A DNAstringSet object or the path to a fasta file
#' @param vcf A vcfR or GRanges object to keep only the windows with the same chromosome name as the vcf/GRanges
#' @param windows_size size of the wondows in base pair. Default is 20000bp
#' @param steps number od base pair for slidding windows. NULL by default gives non-overlapping windows
#' @return A GRanges object with windows
#'
#' @author  ~~Benjamin Laenen~~
#'  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' reference_file <- "my_species_ref.fasta"
#' 
#' #Read the reference
#' REF <- readDNAStringSet(reference_file)
#'
#' #Make windows of 20kb non overlapping
#' windows1 <- windows_from_reference(reference_file)
#' windows2 <- windows_from_reference(REF)
#' identical(windows1, windows2)
#' 
#' #Make windows of 50kb overlapping by 10kb
#' windows <- windows_from_reference(REF, windows_size = 50000, steps = 10000)
#' 
#' 
#' @export
windows_from_reference <- function(reference, vcf = NULL, windows_size = 20000, steps = NULL){
	if(is.null(steps)) steps <- windows_size
	#read the reference, for speeding the function
	#we only select the chr present in the vcf file
	if(class(reference) == "DNAStringSet"){
		REF <- reference
	}else{
		REF <- readDNAStringSet(reference)
	}
	names(REF)  <- gsub(" ", "", names(REF) )
	if(!is.null(vcf) & class(vcf) =="vcfR") REF <- REF[names(REF) %in% unique(getCHROM(vcf))]
	if(!is.null(vcf) & class(vcf) =="GRanges") REF <- REF[names(REF) %in% unique(vcf@seqnames)]

	#creating a Grange object with the genome
	gr <- GRanges(
	    seqnames=Rle(names(REF)),
	    ranges=IRanges(rep(1, length(REF)), end=width(REF)),
	)
	windows_genome <- slidingWindows(gr, windows_size, step = steps)
	#windows_genome <- lapply(gr, GenomicRanges::slidingWindows, windows_size, step = steps)
	return(windows_genome)
}


#' Expand all ranges from a Granges
#' 
#' Function that will create a new Granges object 
#' with all the position from the ranges.
#' 
#' It can be slow for large Grange and should be used
#' if very nescesary. Otherwise more efficent for
#' Genomic Ranges should be used
#' 
#' 
#' @inheritParams remove_overlapping_range
#' @return A GRanges will all range exanped with width of 1
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' 
#' Grange <- GRanges(seqnames="chr1", ranges= c(1:10, 1:4, 2:5, 1:2, 100:1000, 120:140)) 
#' all_sites_Grange <- expand_all_pos_Granges(Grange)
#' #check that all ranges have been expanded
#' # if there was duplicated or overlapping in Grange
#' there will be duplicated in the results
#' identical(length(all_sites_Grange), sum(width(Grange)))
#' sum(duplicated(all_sites_Grange))
#' 
#' #first merge range in the Grange to
#' #avoid duplicates
#' Grange <- reduce(Grange)
#' all_sites_Grange <- expand_all_pos_Granges(Grange)
#' identical(length(all_sites_Grange), sum(width(Grange)))
#' sum(duplicated(all_sites_Grange))
#' 
#' @export
expand_all_pos_Granges <- function(Grange) {
	#Function that will create a new Granges
	#object with all the position from the ranges
	#it can be slow for large Grange
	S <- start(Grange)
	E <- end(Grange)
	POS <- unlist(sapply(1:length(S), function(i) S[i]:E[i]))
	chr <- rep(as.character(seqnames(Grange)), E-S+1)
	res <- GRanges(paste(chr, paste(POS, POS, sep  ="-"), sep = ":"))
	return(res)
}


#' Determine the file separator of a table
#' 
#' The function read the first line of the file
#' and check witch separator is used among common
#' separators : , : ; tab space
#' 
#' @param file character : path to a file
#' @return The separator as a character
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' determine_file_separator("mydata.csv")
#' 
#' @export
determine_file_separator <- function(file){
	#determining separator from the first line of a file
	header <- readLines(file, n=1)
	separators <- c(",",":",";","\t"," ")
	sep <- which(sapply(separators, function(x) grepl(x, header)))
	return(separators[sep])
}

#' A function to check if element of a vector are upper of lower case
#' 
#' A function to check if elemeent are upper of lower case. If the string
#' vector is to long a quick check can be used on a random sample of elements to
#' check there case. Note that mixed case will always return FALSE (e.g. "Ab"),
#' 
#' 
#' @param string A vector of string to check
#' @param Sample should the string vector be subsample to do a quick test assuming that the case (upper or lower) is the same for all element of the vector. Default is NULL, all element are checked.
#' @param isupper Test if string is upper is TRUE, lower if FALSE
#' @return A logical vector of the size of string or size of subsample
#'
#' @author  ~~Benjamin Laenen~~
#'  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' check_case(LETTERS[1:10])
#' check_case(c(LETTERS[1:10], letters[1:10]))
#' check_case(c(LETTERS[1:10], letters[1:10]), isupper = FALSE)
#' check_case(c(LETTERS[1:10], letters[1:10]), Sample = 5, isupper = FALSE)
#' @export
check_case <- function(string, Sample = NULL, isupper = TRUE){
	if(!is.null(Sample)){
		if(length(string) < Sample) Sample = length(string)
		string <- sample(string, Sample)
	}
	if(isupper){
		res <- grepl("^[[:upper:]]+$", string)
	}else{
		res <- grepl("^[[:lower:]]+$", string)

	}
	#determining separator from the first line of a file
	return(res)
}

#' Try mclapply
#' 
#' This function will first try mclapply using multicore. Sometimes it fails
#' without a warnings and return NULL. For the failed run, it will try to run it
#' one signle core using lapply and report a warnings. If lapply fail the code
#' break and an error is thrown.
#
#' 
#' 
#' @param X a vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by as.list.
#' @param FUN the function to be applied to (mclapply) each element of X parallel
#' @param ... see mclapply for details
#' @return  a list of the same length as X and named by X.
#'
#' @author  ~~Benjamin Laenen~~
#'  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' 
#' 
#' @export
mclapply_try <- function(X, FUN, use_BiocParallel = FALSE, no_cores = detectCores(), ...){
	if(!use_BiocParallel){	
		tryout <- mclapply(X, FUN, mc.cores = no_cores ,...)
		#retry if some fails
		failed <- sapply(tryout, is.null)
		if(any(failed)){
			message(sprintf("Run %s failed, trying using single core", paste(which(failed), collapse = " ")))
			tryout[failed] <- lapply(X[failed], FUN)
		} 	
	}else{
		#param <- SnowParam(type = "SOCK", workers = no_cores, ...)
		param <- MulticoreParam(workers = no_cores, ...)
		# bod <- body(FUN)
		# if (trimws(as.character(bod[[1]])) == "{"){
		#     body(FUN)[[2]] <- quote(suppressPackageStartupMessages(library(filterVCF)))
		#     if (length(bod) > 1) body(FUN)[3:(1+length(bod))] <- bod[-1]
		# } else {
		#     body(FUN)[[1]] <- as.name('{')
		#     body(FUN)[[2]] <- quote(suppressPackageStartupMessages(library(filterVCF)))
		#     body(FUN)[[3]] <- bod
		# }
		# fun <- function(X, FUN, param){
		# 	suppressPackageStartupMessages(library(filterVCF))
		# 	bplapply(X, FUN, BPPARAM = param)
		# }
		# tryout <- fun(X, FUN, param = param)
		tryout <- bplapply(X, FUN, BPPARAM = param)
	}
	return(tryout)
}

#' intersect and subtract Granges
#' 
#' This function is inspired from bedtools intersect and bedtools subtract and
#' try to provide a simple way to intersect ranges. vcfR object can be used as
#' either Granges1 or Granges2 but a vcfR object is only return when the vcf is
#' specified for Grange1. The metadata store in Granges1 will be lost if the
#' option merge_interval is TRUE as there is no rules to how to merge the
#' information.
#' 
#' This function is a simple implementation, for more complicated intersection 
#' and customization see the package HelloRanges.
#' 
#' NOTE!:The function ignore strand for the moment
#' 
#' @param Grange1 GRanges, vcfR or path to a BEDfile
#' @param Grange2 GRanges, vcfR or path to a BEDfile
#' @param invert = FALSE If TRUE subtract Grange1 from Grange2
#' @param merge_interval = TRUE Should the interval be merged into larger, note that metadata are lost if TRUE
#' @return intersected ranges
#' \item{ }{if Grange1 is vcfR returns a vcfR} \item{ }{if Grange1 is GRanges returns a GRanges}  \item{ }{if if Grange1 is character returns a GRanges}
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' Grange1 <- GRanges(seqnames="chr1", ranges= c(1:10, 1:4, 2:5, 1:2, 100:1000, 120:140)) 
#' Grange2 <- GRanges(seqnames="chr1", ranges= c(1:10, 60:70, 100:110, 120:140)) 
#'
#' intersect_Granges(Grange1, Grange2)
#' intersect_Granges(Grange1, Grange2, invert = TRUE)
#' intersect_Granges(Grange1, Grange2, merge_interval = FALSE)
#' 
#' @export
intersect_Granges <- function(Grange1, Grange2, invert = FALSE, merge_interval = TRUE){
	# !!!!!!! IMPORTANT
	# UPDATE : now I found the solution in the package HelloRanges
	# I dont not load it as I need only one function but it might
	# be useful in the future. 
	# The function NOW behaves exactly as bedtools intersect and return
	# a vcf file if a vcf is provided as first arguments
	
	# The intersection used to only retain windows that intersect and
	# do not break those windows
	# it is fine if a vcf is Grange1
	# but not when we want to get intersection of 2 windows
	# Solution 1: use expand Grange -> do the intersection -> reduce
	# Solution 2: find the function in Genomicranges that can do that!
    # pairs <- findOverlapPairs(Grange1, Grange2, ignore.strand = TRUE)
    # ans <- pintersect(pairs, ignore.strand = TRUE)


	if(class(Grange1) == "vcfR"){
		vcf <- Grange1
		isVCF <- TRUE
		Grange1 <- vcf2Grange(Grange1)
	}else{
		if(class(Grange1) == "character"){
			if(isTRUE(file.exists(Grange1))) Grange1 <- bed2Grange(Grange1)
		}
		isVCF = FALSE
	}
	if(class(Grange2) == "vcfR"){
		Grange2 <- vcf2Grange(Grange2)
	}else{
		if(class(Grange2) == "character"){
			if(isTRUE(file.exists(Grange2))) Grange2 <- bed2Grange(Grange2)
		}
	}

	if(isTRUE(isVCF)){
		if(!invert){
			filters <- !is.na(GenomicRanges::findOverlaps(Grange1, Grange2 , select ="first", ignore.strand=TRUE))
		}else{
			filters <- is.na(GenomicRanges::findOverlaps(Grange1, Grange2 , select ="first", ignore.strand=TRUE))
		}
		return(vcf[filters])
	}else{
		if(!invert){
		    pairs <- findOverlapPairs(Grange1, Grange2, ignore.strand = TRUE)
	    	ans <- pintersect(pairs, ignore.strand = TRUE)
		}else{
		    hits <- findOverlaps(Grange1, Grange2, ignore.strand = TRUE)
		    toSubtract <- reduce(extractList(Grange2, as(hits, "List")),ignore.strand = TRUE)
		    ans <- psetdiff(Grange1, toSubtract, ignore.strand = TRUE)
		    ans <- subset(unlist(ans), width > 0L)
		    ans
		}
		elementMetadata(ans) <- elementMetadata(ans)[,colnames(elementMetadata(ans)) != "hit"]
		if(merge_interval){
			return(reduce(ans))
		}else{
			return(ans)
		}		
	}
}

#' Get a genomic range with REF and ALT bases
#' 
#' The functions extract the REF and ALt column present in the original vcf
#' file and report the base. It ignores the invariant (ALT = "*") but it
#' includes the fixed sites.
#' 
#' 
#' @param filterVCF a filterVCF object
#' @return A genomicRange object with position from the vcf and REF and ALT
#' base as metadata
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' ref_alt <- get_alternate_GRange(filterVCF)
#' summary(elementMetadata(ref_alt))
#' @export
get_alternate_GRange <- function(filterVCF){
  alt <- vcf2Grange(vcfRaw(filterVCF))
  elementMetadata(alt)["REF"] <- getREF(vcfRaw(filterVCF))
  elementMetadata(alt)["ALT"] <- getALT(vcfRaw(filterVCF))
  alt <- intersect_Granges(alt, c(vcf2Grange(vcf(filterVCF)), fixed(filterVCF)), merge_interval = FALSE)
  elementMetadata(alt) <- elementMetadata(alt)[,1:2]
  return(sort(alt))
}



#' Initialize the option to run filterVCF.R
#' 
#' Should not be used by end user.
#' 
#' Set all the default of opt when initializing 
#' an empty filterVCF object
#' 
#' @return A list of options
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' opt <- parse_args(OptionParser(option_list=option_list))
#' 
#' @export
#' @import optparse
initialise_option <- function(){
	option_list = list(
		make_option(c("-I", "--vcf_file"), action="store", default=NA, type='character', metavar = "file",
		            help=paste(strwrap("Input VCFfile or VCFFile.gz", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),
		
		make_option(c("-R", "--reference"), action="store", default=NA, type='character',
		            help=paste(strwrap("Reference fasta file", width = 80,initial = "\n\t\t", prefix="\n\t\t"), collapse =" ")),
		
		make_option(c("-O", "--output_file"), action="store", default=NA, type='character', metavar = "file",
		            help=paste(strwrap("output name for the filtered VCF.gz. NOTE that if not given output is named automatically based on the input", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-B", "--bed_file"), action="store", default=NA, type='character', metavar = "BED file(s)",
		            help=paste(strwrap("One or multiple bedfile(s) in a list to use for filtering the VCFfile.     [E.g  -B file1.bed,file2.bed,file3.bed]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-r", "--repeats"), action="store", default=NA, type='character', metavar = "BED file",
		            help=paste(strwrap("A bed file with repeats e.g. from RepeatMasker.", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-w", "--filter_repeats_by_windows"), action="store_true", default=FALSE, type='logical', 
		            help=paste(strwrap("Apply a filter to remove 20kb windows that have more than 50% of the bases  covered by repeats defined in the repeat bedfile.  (In practise it can be extended to any bedfile provided in --repeats option). [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-D", "--filter_depth"), action="store", default=NA, type='character', metavar = "min_DP,MAX_DP",
		            help=paste(strwrap("Change the genotypes that don't pass the min and max depth to missing. At least the min depth should be present.     e.g. -D 10,200 or -D 10 [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-g", "--filter_genotype_quality"), action="store", default=NA, type='character', metavar = "GQ,RGQ",
		            help=paste(strwrap("Change the genotypes that don't pass the min GQ/RGQ to missing. If the option is set to auto then the 5% lower quartile is used to determined the threshold. RGQ value are present only in invariable sites. There is no clear recommendation in GATK to filter accroding to GQ or RGQ but to plot and try and use judgment.     e.g. --filter_genotype_quality 5,10 or --filter_genotype_quality 'auto'  [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-G", "--filter_GATK_info"), action="store", default=NULL, type='character', metavar = "GATK like querry",
		            help=paste(strwrap("Apply a filter based on INFO field in the VCFfile populated by GATK. The filters can be applied to any field present in info but recommended are QD, SOR, MQ, MQRankSum, FS, ReadPosRankSum and InbreedingCoeff. Use -G recommended to use GATK recommended threshold for this filter. e.g -G  FS > 60.0 || MQ < 50.0  . Default values are slightly modified from GATK recommendation  https://software.broadinstitute.org/gatk/documentation/article.php?id=6925                           [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-f", "--filter_fix_het"), action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("Apply a filter to remove fixed heterozygous sites in all samples. NOTE : Sites with fixed heterozygous and missing are also removed [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("--filter_all_het"), action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("Apply a filter to remove all heterozygous sites. NOTE that is risky and you have to know what you are doing. This suitable for haploid or hilghy selfing for example. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-m", "--missing"), action="store", default=NA, type='double', metavar = "Proportion",
		            help=paste(strwrap("Apply a filter to remove sites with n proportion of missing data. e.g. -m 0.2 [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-M", "--missing_per_individual"), action="store", default=NA, type='double', metavar = "Proportion",
		            help=paste(strwrap("Apply a filter to remove samples with n proportion of missing data. e.g. --missing_per_individual 0.2 [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-i", "--filter_indel"), action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("Apply a filter to remove indels [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-d", "--filter_high_DP_standardized"), action="store", default=NA, type='character', metavar = "list",
		            help=paste(strwrap("Apply a filter to remove windows that have high standardized (divided by median depth) depth and high depth variance across samples.  Recommended for Arabis is to use 1kb non overlapping windows  and a threshold of >2 or 1.8 for the ratio of median depth or >4 for the 95% quartile of depth. -d (threshold_median1,..,threshold_median2),threshold_CI,windows_size,sliding,percent_sample_passing e.g. : -d (2 1.2 1.5 2.2),4,1000,1000,0.95 If multiple threshold median are set, the first one is used to filter and bedfiles are produce for the others. \n[default %default]", width = 125,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-b", "--biallelic"), action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("apply a filter to keep only bi-allelic sites. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-F", "--filter_fix_het_contiguous_in_pop"), action="store", default=NA, type='character', metavar = "pattern",
		            help=paste(strwrap("Apply a filter to remove contigous  fixed heterozygous sites in a pop by 10bp windows.  A windows is removed if at least 3 fixed hets are found in 10b windows in the pop, which indicates probable a wrongly mapped reads.  The arguments to provide is a regex pattern to match sample names to define the pop e.g. -F 'Riastan' or -F 'ind1|ind2|ind3' or -F '^sample_32.+_.+_pop_[0-9]' \n[default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-p", "--plot"), action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("Produce plots [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--split", action="store", default=NA, type='double', metavar = "integer",
		            help=paste(strwrap("Split the input vcf into n chunks to use multicore. If NA the file is split into the number of cores. If the gz file is larger than 2Gb split the vcf by 500.000 lines. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--gzip_program", action="store", default="unpigz", type='character', metavar = "program name",
		            help=paste(strwrap("Choose which program is used to zip unzip file. unpigz is faster than gzip or gunzip [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--report", action="store", default="vcf_name.reports", type='character', metavar = "file",
		            help=paste(strwrap("Produce a html reports of the filtering. By default it use the vcf name to name the reports. To change the default use --report reports_name.html. To disable the report use --report none", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--single_chrom", action="store_true", default=NA, type='logical',
		            help=paste(strwrap("Does the VCF have only one chromosome? If active skip the checking (quicker), otherwise the first and last line of the VCF are check to see if chr matches. If there is multiple chromosome sand the file (gz) is smaller than 4Gb (500mb) then the splitting is done internally in R and if the file is larger, the VCF is splitted by chromosome and a batch script is written to run the pipeline on each chromosome. It is recommended to split by chromosome for large file to avoid loading large file in R.", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--sample", action="store", default=NA, type='character', metavar = "list or pattern",
		            help=paste(strwrap("Select a list of samples to keep. Use either a list, a file or a regex pattern as follow : --sample sample1,sample2,...,samplen\n--sample sample_file\n--sample pattern='.+Pop1.+' [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--popfile", action="store", default=NA, type='character', metavar = "list or pattern",
		            help=paste(strwrap("A two column tab separated file with sample ID in column 1 and population assignment in column 2. This is used to plot the SFS per population [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-k", "--keep_sites"), action="store", default=NA, type='character', metavar = "file",
	              help=paste(strwrap("BED file with sites to intersect with the vcf.", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		#ad option to selct chr
		make_option("--chr", action="store", default=NA, type='character', metavar = "chromosome name(s)",
		            help=paste(strwrap("EXPERIMENTAL! Select one chromosome or more chromosomes. If only one chromosome is selected, the option single chrom is turn on and the pipeline is faster. e.g --chr chr2 or --chr chr2,ch3,scaffold_300 [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		#add option to have a pop file and analysis is done by pop!
		make_option(c("-c", "--no_cores"), action="store", default=1, type='integer', metavar = "integer",
		            help=paste(strwrap("Number of core [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option(c("-v", "--verbose"), action="store_true", default=TRUE, type='logical',
		            help=paste(strwrap("Print the steps, note that if multicore it is printed to a file called log_filterVCF.R. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),
		
		make_option("--debug", action="store_true", default=NA, type='logical',
		            help=paste(strwrap("Save a .rds object that can be load in R using readRDS. Usefull to recreate the reports using diffent setting for windows size xlab etc. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--merge_only", action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("Used after calling the script with dry_run to run the it with severall batch job. Gather all the rds results, save output and create report. [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n")),

		make_option("--dry_run", action="store_true", default=FALSE, type='logical',
		            help=paste(strwrap("This option allows to split the vcf in chunk and create batch job to run chunks separaterly. A batch job is also created to merge the results [default %default]", width = 120,initial = "\n\t\t", prefix="\n\t\t"), collapse ="\n"))
	)
	return(option_list)

}

#' Merge non overlapping filterVCF object
#' 
#' For computing efficency, a filterVCF can be splitted by
#' chromosome or chunked into lots of piece to run in parallel.
#' This function merge the chunk together. If a chunk missed a slot
#' or is empty and empty slot is merged. The user must ensure that the
#' size of the final vcf, invariant,... equals the sums of the chunked 
#' vcf, invariant. All the slots will be merged except the refrence and 
#' the opt, which are inherited from the first element of the input.
#' 
#' Note that filterVCF object can become quite large in memory when there
#' is many samples sequenced. The option drop_raw_vcf keep only the first
#' sample of the original vcf to free space while maitining the number of
#' sites correct in the original vcf. It means that the original vcf cannot
#' be recreated from a merged object when using drop_raw_vcf=TRUE but the
#' filtered vcf can be outputed with all the fields.
#' 
#' 
#' @param output_multicore list of filterVCF object
#' @param drop_raw_vcf =FALSE Keep only the first sample in the original vcf
#' @return A merged filterVCF object
#'
#' @author  ~~Benjamin Laenen~~
#'
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' 
#' 
#' #out is a list from multiple run of filterVCF.R or multiple chromosomes
#' filtervcf <- merge_output_from_multicore(out, drop_raw_vcf = TRUE)
#' @export
merge_output_from_multicore <- function(output_multicore, drop_raw_vcf = FALSE){
	#check if it exist
	for(i in seq_along(output_multicore)){
		if(class(output_multicore[[i]]) == "list") output_multicore[[i]] <- as.filterVCF(output_multicore[[i]])
	}
	
	#create an empty filterVCF object
	res <- filterVCF()
	reference(res) <- reference(output_multicore[[1]])
	
	if(!drop_raw_vcf){
		message("merging original SNP VCF")
		VCF <- sapply(output_multicore, vcfRaw)
		vcfRaw(res) <- do.call("rbind", VCF)		
	}else{
		message("not merging original SNP VCF but keeping only one sample to be able to create the report")
		VCF <- sapply(output_multicore, function(x){
			vcfRaw(x)@gt <- vcfRaw(x)@gt[,1:2]
			return(vcfRaw(x))
		})
		vcfRaw(res) <- do.call("rbind", VCF)		
	}
	rm(VCF)

	message("merging filtered SNP VCF")
	vcf(res) <- do.call("rbind", sapply(output_multicore, vcf))

	message("merging GRanges invariable")
	invariantRaw(res) <- reduce(do.call("c", sapply(output_multicore, invariantRaw)))

	message("merging GRanges invariable filtered")
	invariant(res) <- reduce(do.call("c", sapply(output_multicore, invariant)))

	message("merging GRanges fixed invariant")
	fixed(res) <- reduce(do.call("c", sapply(output_multicore, fixed)))

	message("merging master filter")
	masterFilter(res)  <- do.call("c", lapply(output_multicore, masterFilter))

	message("merging master filter for invariable")
	masterFilterInv(res) <- do.call("c", lapply(output_multicore, masterFilterInv))

	message("merging single filter")
	filters <- lapply(output_multicore, Filters)
	filters_name <- names(filters[[1]])
	Filters(res) <- lapply(filters_name,function(x) do.call("c", lapply(filters, "[[", x)))
	names(Filters(res)) <- filters_name
	rm(filters)

	message("merging single filter for invariable")
	filters_inv <- lapply(output_multicore, FiltersInv)
	filter_inv_name <- names(filters_inv[[1]])
	FiltersInv(res) <- lapply(filter_inv_name,function(x) do.call("c", lapply(filters_inv, "[[", x)))
	names(FiltersInv(res)) <- filter_inv_name
	rm(filters_inv)

	message("merging GRanges for removed sites")
	removed_sites(res) <- do.call("c", lapply(output_multicore, removed_sites))

	message("merging GRanges for removed sites invariable")
	removed_sites_inv(res) <- do.call("c", lapply(output_multicore, removed_sites_inv))

	if(is.null(BEDFilter(output_multicore[[1]]))){
		BEDFilter(res) <- GRanges()
	}else{
		BEDFilter(res) <- reduce(do.call("c", lapply(output_multicore, BEDFilter)))
	}

	message("merging GRanges for windows with half repeats to remove")
	if(is.null(Filter_WindowsProp(output_multicore[[1]]))){
		Filter_WindowsProp(res) <- GRanges()
	}else{
		Filter_WindowsProp(res) <- reduce(do.call("c", lapply(output_multicore,Filter_WindowsProp)))
	}

	message("merging GRanges for fix heterozygous")
	if(is.null(FixHet_RefPop(output_multicore[[1]]))){
		FixHet_RefPop(res) <- GRanges()
	}else{
		FixHet_RefPop(res) <- reduce(do.call("c", lapply(output_multicore, FixHet_RefPop)))
	}

	message("merging GRanges for filter for repeats by depth")
	filter_DP_repeats <- lapply(output_multicore, Normalized_DP_Filter)
	filter_DP_repeats_name <- names(filter_DP_repeats[[1]])
	filter_DP_repeats_1 <- lapply(filter_DP_repeats, "[[", filter_DP_repeats_name[1])
	filter_DP_repeats_2 <- lapply(filter_DP_repeats, "[[", filter_DP_repeats_name[2])
	filter_DP_repeats_1 <- lapply(seq_along(filter_DP_repeats_1[[1]]),function(x) do.call("c", lapply(filter_DP_repeats_1, "[[", x)))
	filter_DP_repeats_2 <- lapply(seq_along(filter_DP_repeats_2[[1]]),function(x) do.call("c", lapply(filter_DP_repeats_2, "[[", x)))
	filter_DP_repeats <- list()
	filter_DP_repeats[[1]] <-filter_DP_repeats_1
	filter_DP_repeats[[2]] <-filter_DP_repeats_2
	names(filter_DP_repeats) <- filter_DP_repeats_name
	Normalized_DP_Filter(res) <- filter_DP_repeats
	
	res@stats_vcf_filtered <- NA_character_

	GetOpt(res)  <-  GetOpt(output_multicore[[1]])

	rm(output_multicore)
	gc(verbose = FALSE)
	return(res)
}

#' Subsample a filterVCF object
#' 
#' Extract a subsample using a list of individuals and recalculate SNPs, invariant and fixed sites
#' 
#' 
#' @param filtervcf a filterVCF object
#' @param subsample regex pattern to grep the VCF column names. Can be a list of individual separated by "|"
#' @param missing_allowed Percentage of missing data allowed in the subsample
#' @param gt a genotype matrix extracted using extract_gt_ff(), speeds up the function
#' @return a subsampled filterVCF object
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~filterVCF
#' @examples
#' 
#' sub_filterVCF  <- subsample_filterVCF(filterVCF_object, subsample = "ind1|ind3|ind234")
#' 
#' @export

subsample_filterVCF <- function(filtervcf, subsample, missing_allowed = NA, gt = NULL){
	#only work for one chromosome, use selectCHROM then run the subsample
	#subsample = "ITM2"
	all_sample <- colnames(vcf(filtervcf)@gt)[-1]
	sample_list <- sort(grep(subsample, all_sample, value = TRUE))

	#If there is no sample we cannot subsample!
	if(length(sample_list) == 0){
		message(sprintf("Samples %s are not present in the VCF, check that names in the popfile match names in the VCF"))
		return(filterVCF())
	}

	# Take only the subsample from the filtered VCF
	VCF <- vcf(filtervcf)[,c("FORMAT", sample_list)]
	

	if(is.null(gt)){
		gt <- extract_gt_ff(VCF)
	}

	sub_filtervcf <- filterVCF()
	#index of sites that are invariant
	#in the filtered VCF
	#We cannot work on the VCFRaw because sometime we condense it by 
	#only retaining one sample thus we cannot subsample.
	filter_invariant <- create_invariant_filter(VCF, gt=gt)
	#get the index of the VCFRaw that correspond to the invariant
	index_inv_vcfRaw <- which(getPOS(vcfRaw(RES)) %in% getPOS(VCF[filter_invariant]))

	filter_fixed <- create_filter_fixed_sites(VCF, gt=gt)
	#get the index of the VCFRaw that correspond to the invariant
	index_fixed_vcfRaw <- which(getPOS(vcfRaw(RES)) %in% getPOS(VCF[filter_fixed]))

	filter_missing <- create_filter_missing(gt=gt,  allowed_missing_threshold =  if(is.na(missing_allowed)) GetOpt(filtervcf)$missing else missing_allowed)	
	index_missing_vcfRaw <- which(getPOS(vcfRaw(RES)) %in% getPOS(VCF[!filter_invariant & !filter_fixed & filter_missing]))
	index_missing_Newinvariant <- which(getPOS(vcfRaw(RES)) %in% getPOS(VCF[filter_invariant  & !filter_fixed & filter_missing]))
	index_missing_NewFixed <- which(getPOS(vcfRaw(RES)) %in% getPOS(VCF[!filter_invariant  & filter_fixed & filter_missing]))

	vcfRaw(sub_filtervcf) <- vcfRaw(RES)[-unique(c(index_inv_vcfRaw, index_fixed_vcfRaw))]

	#the new filtered VCF
	#minus invairant minus fixed and minus missing the sub sample
	vcf(sub_filtervcf) <- VCF[!filter_invariant & !filter_fixed & !filter_missing,]
	
	#Invariant raw
	# + new invariant
	invariantRaw(sub_filtervcf) <- reduce(c(invariantRaw(RES), vcf2Grange(VCF[filter_invariant])))

	#invariant filtered
	# - missing
	invariant(sub_filtervcf) <- reduce(c(invariant(RES), vcf2Grange(VCF[filter_invariant & !filter_missing])))

	#fixed
	# + fixed - missing
	fixed(sub_filtervcf) <- reduce(c(fixed(RES), vcf2Grange(VCF[filter_fixed & !filter_missing])))
	
	#remove sites from the SNP
	# +missing - invariant - fixed
	removed_sites(sub_filtervcf) <-  reduce(c(removed_sites(RES), vcf2Grange(VCF[!filter_invariant & !filter_fixed & filter_missing])))

	#remove sites from the inv
	# +missing - invariant - fixed
	removed_sites_inv(sub_filtervcf) <-  reduce(c(removed_sites_inv(RES), vcf2Grange(VCF[filter_invariant & filter_fixed & filter_missing])))

	#master filter
	# the size changes to the new length of vcfRaw
	# + missing
	masterFilter(sub_filtervcf) <- masterFilter(RES)
	masterFilter(sub_filtervcf)[index_missing_vcfRaw] <- TRUE
	masterFilter(sub_filtervcf) <- masterFilter(sub_filtervcf)[-unique(c(index_inv_vcfRaw, index_fixed_vcfRaw))]

	#master filter  inv
	# the size changes to the new length of vcfRaw
	# + missing
	# we never use this one as we dont have the
	# original VCF anymore
	masterFilterInv(sub_filtervcf) <- logical(0)

	#filters
	# add to the missing
	# change the size of all filters 
	# to remove the invariant
	Filters(sub_filtervcf) <- Filters(RES)
	Filters(sub_filtervcf)$missing[index_missing_vcfRaw] <- TRUE 
	Filters(sub_filtervcf) <- lapply(seq_along(Filters(sub_filtervcf)), function(i) Filters(sub_filtervcf)[[i]][-unique(c(index_inv_vcfRaw, index_fixed_vcfRaw))])
	names(Filters(sub_filtervcf)) <- names(Filters(RES))


	#filters inv, approximation here
	#I have no way of going back
	#then I cannot add the missing in the pop
	#and the new fixed or new invariant
	#to those filters
	#The difference will be the sites that are only invariant or fixed in the pop
	#Those can be retries in the invariant and fixed slot
	#by comparing be fore and after the subsampling
	# add to the missing
	# change the size of all filters 
	# to remove the invariant
	FiltersInv(sub_filtervcf) <- FiltersInv(RES)

	reference(sub_filtervcf) <- reference(RES)
	GetOpt(sub_filtervcf) <- GetOpt(RES)

	Normalized_DP_Filter(sub_filtervcf) <- Normalized_DP_Filter(RES)
	BEDFilter(sub_filtervcf) <- BEDFilter(RES)
	FixHet_RefPop(sub_filtervcf) <- FixHet_RefPop(RES)
	Filter_WindowsProp(sub_filtervcf) <- Filter_WindowsProp(RES)

	message(sprintf("\nKeeping %s samples:\n%s", length(sample_list), paste(sample_list, collapse = "\n")))
	message(sprintf("\nRemoving %s sites with missing data in the subsample", sum(filter_missing & !filter_invariant & !filter_fixed)))
	message(sprintf("\n%s sites became invariant in the subsample", sum(filter_invariant)))
	message(sprintf("\n%s sites among sites that became invariant are missing in the subsample", sum(filter_invariant & filter_missing)))
	message(sprintf("\n%s sites became fixed in the subsample", sum(filter_fixed)))
	message(sprintf("\n%s sites among sites that became fixed are missing in the subsample", sum(filter_fixed & filter_missing)))
	message(sprintf("\nRemaining %s variant sites in the subsampled vcf", nrow(vcf(sub_filtervcf))))
	return(sub_filtervcf)
}