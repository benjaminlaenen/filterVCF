# ==========================================================================
#                  Vcf utlities
# ==========================================================================


#' Utility to split a VCF file into chunks
#' 
#' This function is used internally in the pipeline filterVCF.R but can be called by user.
#' The function will split any vcf file into equal chunks depending on the
#' nb_files argument. However if the gziped VCF file is larger than 2Gb then
#' chunks of 5.000.000 lines are created by default. This function is intended
#' to be used with very large VCF (20-100Gb) file resulting from a calling of
#' all sites in the genome. Note however that is not compatible with the GVCF
#' format of GATK with ranges to represent invariant chunks.
#' 
#' NOTE only tested on UNIX system
#' 
#' 
#' @param vcf_file path to a vcf file (gziped or not)
#' @param outputdir Output directory to save the splitted vcf
#' @return a character vector with the names of the splitted files.
#' 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' opt = initialise_option()
#' opt$gzip_program <- "gzip"
#' opt$no_cores <- 4
#' split_vcf_multiple_file("my.vcf", "split_out", opt)
#' 
#' @export
split_vcf_multiple_file <- function(vcf_file,  outputdir="./", nb_files=NA, ...){
	#Function to split a vcf file into n parts into the current directory or a new
	#one if specified.
	#It can handle gzip vcf and will output a gzip vcf.
	#Use multicore if specified in opt
	# .libPaths("/proj/uppstore2018024/private/Rpackages/")
	# suppressPackageStartupMessages(suppressMessages(try(library(filterVCF))))
	myDots <- list(...)
	 if (!is.null(myDots$opt)){
	    opt <- myDots$opt
	 }else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	 }
	vcf_size <- file.size(vcf_file)
	if(is.na(nb_files)) nb_files <- opt$no_cores

	if(isGzipped(vcf_file)){
		options(scipen=999)
		require(parallel)
		if(opt$verbose) write(strwrap(sprintf("VCF  %s.vcf is a gzip file using %s to decompress.\n", basename(vcf_file), opt$gzip_program), width=80),stderr())

		outputdir <- gsub(".gz", "", outputdir)
		vcf_file <- basename(vcf_file)
		vcf_file <- gsub(".vcf.gz", "", vcf_file)

		# if(Sys.info()['sysname'] == "Darwin"){

		if(vcf_size > 2e9){
			if(opt$verbose) write(strwrap(sprintf("VCF  %s.vcf is larger than 2Gb and will be splitted by 5.000.000 lines to avoid gunzipping twice .\n", vcf_file), width=80),stderr())
			system(sprintf("%s -c  %s.vcf.gz |  split -l %s - %s", opt$gzip_program, vcf_file, 5000000, vcf_file))
		}else{
			if(Sys.info()["sysname"] == "Darwin"){
				nb_of_line <- strsplit(system(sprintf("%s -c %s.vcf.gz | wc -l", opt$gzip_program,  vcf_file), intern=TRUE), " +|\t")[[1]][2]
			}
			if(Sys.info()["sysname"] == "Linux"){
				nb_of_line <- strsplit(system(sprintf("%s -c %s.vcf.gz | wc -l", opt$gzip_program,  vcf_file), intern=TRUE), " +|\t")[[1]]
			}
			lines_parts <- ceiling(as.numeric(nb_of_line) / nb_files)
			if(opt$verbose) write(strwrap(sprintf("VCF  %s.vcf has %s lines and will be splitted into %s parts with each part having ~ %s lines each.\n", vcf_file,nb_of_line, nb_files, lines_parts), width=80),stderr())
			system(sprintf("%s -c  %s.vcf.gz |  split -l %s - %s", opt$gzip_program, vcf_file, lines_parts, vcf_file))
		}
		if(opt$verbose) message("Splitting done...\n", appendLF=TRUE)

		#extracting the header
		system(sprintf("%s -c  %s.vcf.gz | head -n 30000| grep '^#'  > %s_header", opt$gzip_program, vcf_file, vcf_file))

		FL <- list.files(pattern=paste0(vcf_file, "..$"))

		void <- file.rename(FL[1], paste0(outputdir, "/", vcf_file, "_1.vcf"))
		for(i in seq_along(FL[-1])){
			system(sprintf("cat %s_header %s > %s/%s_%s.vcf", vcf_file, FL[i+1], outputdir, vcf_file, as.character(i+1)) )
		}

		void <- file.remove(FL[-1])
		void <- file.remove(paste0(vcf_file, "_header"))
	}else{
		options(scipen=999)
		vcf_file <- basename(vcf_file)
		vcf_file <- gsub(".vcf", "", vcf_file)

		if(Sys.info()['sysname'] == "Darwin"){
			nb_of_line <- strsplit(system(sprintf("wc -l %s.vcf", vcf_file), intern=TRUE), " +|\t")[[1]][2]
			lines_parts <- round(as.numeric(nb_of_line) / nb_files)
			if(opt$verbose) write(strwrap(sprintf("VCF  %s.vcf has %s lines and will be splitted into %s parts with each part having ~ %s lines each.\n", vcf_file,nb_of_line, opt$no_cores-1, lines_parts), width=80),stderr())
			system(sprintf("split -l %s %s.vcf %s", lines_parts, vcf_file, vcf_file))
			if(opt$verbose) message("Splitting done...", appendLF=TRUE)
		}

		if(Sys.info()['sysname'] == "Linux"){
			if(opt$verbose) write(strwrap(sprintf("VCF  %s.vcf will be splitted into %s parts.\n", vcf_file, opt$no_cores-1), width=80),stderr())
			system(sprintf("split -n %s %s.vcf %s", nb_files, vcf_file, vcf_file))
			if(opt$verbose) message("Splitting done...", appendLF=TRUE)
		}
		#Can cause problem if the header is larger than 100000
		system(sprintf("head -n 100000 %s.vcf | grep '^#' > %s_header", vcf_file, vcf_file))
		FL <- list.files(pattern=paste0(vcf_file, "..$"))
		#change name of the first file that contains the header
		void <- file.rename(FL[1], paste0(outputdir, vcf_file, "_1.vcf"))

		#add the header to all subsquent chunk
		for(i in seq_along(FL[-1])){
			system(sprintf("cat %s_header %s > %s%s_%s.vcf",vcf_file, FL[i+1], outputdir, vcf_file, as.character(i+1)) )
		}
		void <- file.remove(FL[-1])
		void <- file.remove(paste0(vcf_file, "_header"))
	}
	#grep the file name of the splitted vcfs
	names_of_split_file <- list.files(path=getwd(), pattern=paste0(vcf_file, "_[0-9]{+}\\.vcf$"), full.names = TRUE, include.dirs = TRUE, recursive = TRUE)

	if(opt$verbose) message(paste("Gzipping", names_of_split_file, "...\n"), appendLF=TRUE)

	if(opt$no_cores > 1){
		cluster <- parallel::makeCluster(opt$no_cores, type="PSOCK")
		clusterExport(cl=cluster, list("names_of_split_file", "escape_dropbox_pro_path"), envir=environment())
		out <- parLapply(cluster, names_of_split_file, function(x) system(sprintf("gzip -f --fast %s", escape_dropbox_pro_path(x))))
		stopCluster(cluster)
	}else{
		void <- lapply(names_of_split_file, gzip, overwrite=TRUE)
		}
	if(opt$verbose) message("Gzipping done...\n", appendLF=TRUE)

	#reorder alpha numerically
	names_of_split_file <- names_of_split_file[order(as.numeric(gsub(".+_([0-9]{+}).vcf", "\\1", names_of_split_file)))]
	names_of_split_file <- paste0(names_of_split_file, ".gz")
	return(names_of_split_file)
}



#' Check the number of chromosome in a VCF
#' 
#' This function only works on UNIX system and call tail to check if the chrmomosome name of the first line of the VCF differ from the last line. 
#' 
#' Note that it does not check the entire file.
#' 
#' @param vcf_file a path to a VCF file gzip or not
#' @return TRUE if only one chromosome is present
#' @note  ~~further notes~~ 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' @export
#' @importFrom R.utils isGzipped gzip
check_nb_chrom <- function(vcf_file, ...){
	myDots <- list(...)
	 if (!is.null(myDots$opt)){
	    opt <- myDots$opt
	 }else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	 }
	if(file.size(vcf_file) > 5e8){
		message("The VCFfile(gz) is over 500mb, there will be no checking for multiple chromosomes.\nUse option --single_chrom to specify that the VCF has only one chrom. Otherwise the VCF is splited into chrom and a batch script is written to run the pipeline separately.", appendLF=TRUE)
	}else{
		if(opt$verbose) message("Checking if the VCFfile has multiple chromosomes.", appendLF=TRUE)
		last <- system(sprintf("%s -c %s | tail -n 1",opt$gzip_program, escape_dropbox_pro_path(vcf_file)), intern =TRUE)
		if(isGzipped(vcf_file)){
			zz <- gzfile(vcf_file,'rt')
		}else{
			zz <- vcf_file
		}
		#read first and last line
		first_chr <- scan(zz, what="character", skip = 1, n=1, comment.char="#", quiet = TRUE)
		last_chr <- unlist(strsplit(last, split = "\t"))[1]
		if(first_chr == last_chr){
			return(TRUE)
		}else{
			message(sprintf("At least two chromosomes are different in VCF %s %s", first_chr, last_chr))
			close.connection(file(vcf_file))
			return(FALSE)
		}
	}
}



#' Utility to split a VCF file by chromosome
#' 
#' This function read the VCF file using a quick read.line function and save
#' splitted vcf file by chromosome on disk. It calls system to gunzip/unpigz and
#' is only tested on UNIX system. It is used internally in the pipeline
#' filterVCF.R. If you encounter problem, try to gunzip the vcffile first and then split.
#' 
#' 
#' @param vcf_file path to a vcf file (gziped or not)
#' @param outputdir Output directory to save the splitted vcf
#' @return a character vector with the names of the splitted files.
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' opt = initialise_option()
#' opt$gzip_program <- "gzip"
#' split_vcf_by_chrom("my.vcf", "split_out", opt)
#' 
#' 
#' @export
split_vcf_by_chrom <- function(vcf_file, outputdir, ...){
	my.read.lines2 <- function(fname) {
		s = file.info( fname )$size
		buf = readChar( fname, s, useBytes=T)
		strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
	}

	split_vcf <- function(vcf_file){
		VCF <- my.read.lines2(vcf_file)
		#VCF <- readLines(basename(vcf_file))
		index_header <- grepl("^#", VCF)
		header <- VCF[index_header]
		VCF <- VCF[!index_header]

		chr <- gsub("^([A-Z]*[a-z]*[0-9]+)\t[0-9]*\t.+$", "\\1", VCF)

		vcf_by_chr <- split(VCF, chr)
		names(vcf_by_chr) <- levels(as.factor(chr))

		vcf_by_chr <- lapply(vcf_by_chr, function(x) c(header, x))
		return(vcf_by_chr)
	}

	vcf_copy <- paste0(outputdir, "/", gsub(".gz", "", basename(vcf_file)))

	if(isGzipped(vcf_file)){
		system(sprintf("%s -f -c %s > %s", opt$gzip_program, escape_dropbox_pro_path(vcf_file),  escape_dropbox_pro_path(vcf_copy)))
	}else{
		void <- file.symlink(vcf_file, vcf_copy)
	}

	vcf_copy_by_chr <- split_vcf(vcf_copy)
	#save splitted vcf to file
	for(i in names(vcf_copy_by_chr)){
		writeLines(vcf_copy_by_chr[[i]], con = paste0(gsub(".vcf", "", vcf_copy), "_", i , ".vcf"))
	}
	file.remove(vcf_copy)
	lapply(paste0(gsub(".vcf", "", vcf_copy), "_", names(vcf_copy_by_chr) , ".vcf"), gzip, overwrite=TRUE)
	return(paste0(gsub(".vcf", "", vcf_copy), "_", names(vcf_copy_by_chr) , ".vcf.gz"))
}



#' Split a vcfR by GRange
#' 
#' Use a Grange object as a factor to split a vcfR object. Each ranges will be intersect with the vcfR.
#' 
#' 
#' @param vcfR A vcfR object
#' @param Grange A GRange object
#' @return a list of vcfR object the size of length(Grange)
#' 
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~utilities
#' @examples
#' 
#' #Split a vcfR by gene
#' #' gff <- import("my_species.gff3", format = "gff")
#' gff <- gff[values(gff)$type == "gene"]
#' vcf_by_gene <- split_vcfR(vcfR, gff)
#' 
#' @export
split_vcfR <- function(vcfR, Grange){
	Grange1 <- vcf2Grange(vcfR)
	split_factor <- findOverlaps(Grange1, Grange, select ="first", ignore.strand=TRUE)
	splited_gt <- split(vcfR@gt, as.factor(split_factor))
	splited_fix <- split(vcfR@fix, as.factor(split_factor))
	NCOL_gt <- ncol(vcfR@gt)
	NCOL_fix <- ncol(vcfR@fix)
	DIMNAME_gt <- dimnames(vcfR@gt)
	DIMNAME_fix <- dimnames(vcfR@fix)
	splitted_vcf <- lapply(seq_along(splited_gt), function(i) new("vcfR", meta = vcfR@meta, fix = matrix(splited_fix[[i]], ncol = NCOL_fix, byrow = FALSE, dimnames = DIMNAME_fix), gt = matrix(splited_gt[[i]], ncol = NCOL_gt, byrow = FALSE, dimnames = DIMNAME_gt)))
	names(splitted_vcf) <- names(splited_fix)
	
	#add the intervals of the Granges
	#that have no sites in the vcf.
	#The length of the results should be the length 
	#of the Grange
	missing_Grange <- seq_along(Grange)[!seq_along(Grange) %in% unique(split_factor)]
	empty_vcfR_list <- lapply(seq_along(missing_Grange), function(x) new("vcfR"))
	names(empty_vcfR_list) <- missing_Grange

	all_splitted_vcf <- c(splitted_vcf, empty_vcfR_list)
	all_splitted_vcf <- all_splitted_vcf[order(as.numeric(names(all_splitted_vcf)))]

	return(all_splitted_vcf)
}



