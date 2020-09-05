# ==========================================================================
#                    Parser arguments
# ==========================================================================


#' parse option for GATK filter
#' 
#' Intented to be use with the pileline filterVCF.R
#' #' 
#' ~~ If necessary, more details than the description above ~~
#' 
#' @param filter_GATK_info Character string that follow GATK recomentdation
#' @return list of GATK filter to apply
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$filter_GATK_info <- "QD <  5.0  || QD >   6 || FS > 60.0 || MQ < 50.0 || MQ > 65.0 || MQRankSum < -5 || ReadPosRankSum < -4.0 || SOR > 3 || InbreedingCoeff < -0.2"
#' parse_GATK_filter_option(opt$filter_GATK_info)
#' 
#' 
#' @export
parse_GATK_filter_option <- function(filter_GATK_info){
	#opt$filter_GATK_info <- "QD <  5.0  || QD >   6 || FS > 60.0 || MQ < 50.0 || MQ > 65.0 || MQRankSum < -5 || ReadPosRankSum < -4.0 || SOR > 3 || InbreedingCoeff < -0.2"
	filter_GATK_info <- unlist(strsplit(filter_GATK_info, " \\|\\| "))
	filter_GATK_info <- data.frame(matrix(unlist(strsplit(filter_GATK_info, " +")), ncol =3, byrow = TRUE))
	filter_GATK_info$rule <- paste(filter_GATK_info[,2], filter_GATK_info[,3])
	filter_GATK_info$GATK_filter_name <- filter_GATK_info[,1]

	filter_list <- list()
	for(i in filter_GATK_info$GATK_filter_name){
		filter_list[[i]] <- filter_GATK_info$rule[filter_GATK_info$GATK_filter_name == i]
	}
	return(filter_list)
}



#' parse options for filter using normalized depth
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_filter_high_DP_standardized Character string, see help of filterVCF.R pipeline for more info
#' @return list of of filter threshold and windows size to be used during the filtering
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$filter_high_DP_standardized <- "(2 1.2 1.5 2.2),4,1000,1000"
#' parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)
#' 
#' 
#' @export
parse_filter_high_DP_standardized <- function(opt_filter_high_DP_standardized){
	if(isTRUE(all(!is.na(opt_filter_high_DP_standardized)))){
		if(is.list(opt_filter_high_DP_standardized)){
			return(opt_filter_high_DP_standardized)
		}else{
			#opt_filter_high_DP_standardized <- "(2 1.2 1.5 2.2),4,1000,1000"
			opt_filter_high_DP_standardized <- as.list(strsplit(opt_filter_high_DP_standardized, ",")[[1]])
			res <- list()
			res$threshold <- suppressWarnings(eval(parse(text=paste0("c", gsub(" +", ",", opt_filter_high_DP_standardized[[1]]) ))))
			res$threshold_CI <- as.numeric(opt_filter_high_DP_standardized[[2]])
			res$windows_size <- as.numeric(opt_filter_high_DP_standardized[[3]])
			res$slidding <- as.numeric(opt_filter_high_DP_standardized[[4]])
			if(isTRUE(length(opt_filter_high_DP_standardized) == 5)){
				res$percent_passing_filter <- as.numeric(opt_filter_high_DP_standardized[[5]])
				if(log(res$percent_passing_filter) > 0) {
					message(sprintf("The percent of sample passing the high DP normalized depth filter should be between 0 and 1. Now it is %s", res$percent_passing_filter))
					message("setting it to 95%")
					res$percent_passing_filter <- 0.95
				}
			}else{
				res$percent_passing_filter <- 0.95
			}

			return(res)
		}
	}else{
		return(NA)
	}
}



#' parse options for filter by bed files
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_bed_file Character string, see help of filterVCF.R pipeline for more info
#' @return list of of filter threshold and windows size to be used during the filtering
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$bed_file <- "file1.bed,file2.bed,file3.bed"
#' parse_bed_file_options(opt$bed_file)
#' 
#' @export
parse_bed_file_options <- function(opt_bed_file){
	#opt_bed_file <- "file1.bed,file2.bed,file3.bed"
	opt_bed_file <- unlist(strsplit(opt_bed_file, " ?, ?"))
	#test if bedfile have a path or not
	for(i in seq_along(opt_bed_file)){
		opt_bed_file[i] <- normalizePath(opt_bed_file[i])
	}
	return(opt_bed_file)
}



#' parse options for filter individual depth
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_filter_depth Character string, see help of filterVCF.R pipeline for more info
#' @return list of min and max depth threshold
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$filter_depth <- "10,200"
#' parse_filter_depth(opt$filter_depth)
#' 
#' @export
parse_filter_depth <- function(opt_filter_depth){
	if(!is.na(opt_filter_depth)){
	#opt_filter_depth <- "10,200"
	opt_filter_depth <- as.list(strsplit(opt_filter_depth, " ?, ?")[[1]])
	if(length(opt_filter_depth) == 1){
		opt_filter_depth$max_DP <- 1e6
	}else{
		opt_filter_depth$max_DP <- as.numeric(opt_filter_depth[2])
	}
	opt_filter_depth$min_DP <- as.numeric(opt_filter_depth[1])
	return(opt_filter_depth)
	}else{
		return(NA)
	}
}



#' parse options for filter by genotype quality
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_filter_depth Character string, see help of filterVCF.R pipeline for more info
#' @return list of threshold for filter on GQ and RGQ (invariant)
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$filter_GQ <- "5,10"
#' #OR
#' opt$filter_GQ <- "auto"
#' parse_filter_GQ(opt$filter_GQ)
#' 
#' @export
parse_filter_GQ <- function(opt_filter_GQ){
	if(!is.na(opt_filter_GQ)){
	# opt_filter_GQ <- "5,10"
	# opt_filter_GQ <- "auto"
	if(opt_filter_GQ == "auto"){
		opt_filter_GQ <- list(GQ = "auto", RGQ = "auto")
	}else{
		opt_filter_GQ <- as.list(as.numeric(strsplit(opt_filter_GQ, " ?, ?")[[1]]))
		names(opt_filter_GQ) <- c("GQ", "RGQ")
	}	
	return(opt_filter_GQ)
	}else{
		return(NA)
	}
}

#' parse options for selecting samples
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_filter_depth Character string, see help of filterVCF.R pipeline for more info
#' @return regex expression to grep selected samples
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$sample <- "pattern='pop1_'"
#' parse_opt_sample(opt$opt$sample)
#' 
#' @export
parse_opt_sample <- function(opt_sample, ...){
	myDots <- list(...)
	 if (!is.null(myDots$opt)){
	    opt <- myDots$opt
	 }else{
	 	opt <- initialise_option()
	 }
	if(!is.na(opt_sample)){
		if(grepl("pattern", opt_sample)){
			pattern <- gsub("pattern *= *", "", opt_sample)
			#grep the sample name from the header
			#this will work better than reading x line
			#if there is more than x scaffold
			sample <- colnames(read.vcfR(opt$vcf_file, nrow=1, verbose=FALSE)@gt)[-1]
			res <- grep(pattern, sample, value = TRUE)
		}
		if(file.exists(opt_sample)){
			res <- unlist(readLines(opt_sample))
		}
		if(!file.exists(opt_sample) && grepl(",", opt_sample)){
			res <- unlist(strsplit(opt_sample, " *, *"))
		}
		if(!file.exists(opt_sample) && grepl("|", opt_sample)){
			res <- unlist(strsplit(opt_sample, " *\\| *"))
		}
		if(isTRUE(opt$verbose)) message(sprintf("The following samples are kept : %s", paste(res, collapse = "\n\t")), appendLF=TRUE)
		return(res)
	}else{
		return(".")
	}
}



#' parse options for naming the report
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_report Character string, see help of filterVCF.R pipeline for more info
#' @param vcf_file Path to a vcf file, see help of filterVCF.R pipeline for more info
#' @return Character string with report name
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$report <- "vcf_name.reports.html"
#' parse_opt_report(opt$report, vcf_file = "species1.vcf" )
#' 
#' @export
parse_opt_report <- function(opt_report, vcf_file){
	#opt_report="vcf_name.reports.html"
	if(grepl("vcf_name", opt_report)){
		opt_report <- gsub("vcf_name", gsub(".vcf|.vcf.gz","", basename(vcf_file)), opt_report)
	}
	if(opt_report == "none"){
		opt_report <- NA
	}
	return(opt_report)
}



#' parse options for selecting chromosome
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_chr Character string, see help of filterVCF.R pipeline for more info
#' @return Character string with chromosome name
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$chr <- "chr1,chr2, chr3"
#' parse_opt_chr(opt$chr)
#' 
#' @export
parse_opt_chr <- function(opt_chr){
	opt_chr <- unlist(strsplit(opt_chr, " ?, ?")[[1]])
	return(opt_chr)
}

#' parse options for file with population assignment
#' 
#' Intented to be use with the pileline filterVCF.R
#' 
#' @param opt_popfile path to a tab-separated popfile, see help of filterVCF.R pipeline for more info
#' @return dataframe with sample names and population
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords parser
#' @examples
#' 
#' Parse the option from a call to 
#' option_list <- initialise_option()
#' opt <- parse_args(OptionParser(option_list=option_list))
#' opt$popfile <- "popfile.tsv"
#' parse_opt_popfile(opt$popfile)
#' 
#' @export
parse_opt_popfile <- function(opt_popfile){
	if(file.exists(opt_popfile)){
		pop <- read.table(opt_popfile, sep = "\t", stringsAsFactors = FALSE)
		return(pop)
	}else{
		message(sprintf("File %s does not exist, no populations defined", opt_popfile))
		return(NA)
	}
}
