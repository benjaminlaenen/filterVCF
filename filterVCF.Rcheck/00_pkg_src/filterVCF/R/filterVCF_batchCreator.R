# ==========================================================================
#            Batch creator
# ==========================================================================


#' Write a batch script to run filterVCF slurm job
#' 
#' This function should not be used directly but is part of the filterVCF.R
#' pipeline
#' 
#' 
#' @param vcf_file path to the vcf
#' @param outputdir output directory to save the results
#' @param ... list of option produced by initialise_option()
#' @return Write batch job files
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~batchCreator
#' @examples
#' 
#' 
#' @export
write_batch_script_slurm <- function(vcf_file, outputdir,...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}

	bin_bash <- sprintf("#!/bin/bash -l
#SBATCH -A snic2017-7-175
#SBATCH -J filter_vcf
#SBATCH -p core -n %s
#SBATCH -t %s

cd %s
", opt$no_cores,  "60:00:00" , outputdir)
	#cmd <- "Rscript test.R --missing 0.2 -I a --gzip_program gunzip"
	cmd <- get_call()
    # cmd <- call
	cmd <- gsub("--vcf_file '.+.vcf.?.?.?' +-|-I '.+.vcf.?.?.?' +-", paste0("--vcf_file " ,basename(vcf_file), " -"), cmd)
	cmd_keep_rds_only <- ""
	if(opt$dry_run){
		cmd <- gsub("--dry_run", "--split 1", cmd)
		vcf_name <- gsub(".vcf.gz|.vcf", "", basename(vcf_file))
		cmd_keep_rds_only <- sprintf("\n# mv %s.rds ../%s.rds\n\n# cd ..\n\n# rm -r %s", vcf_name, vcf_name, vcf_name)
	}
	message("setting --no_cores to 1, multithreading is not yet optimal, and sometimes create memmory allocation problem. Change it in the batch script to try multicore")
	cmd <- gsub("--no_cores .+ -", "--no_cores '1' -", cmd)
	bash_job_file <- sprintf("Run_filterVCF_%s.sh", gsub(".vcf.gz|.vcf", "", basename(vcf_file)))
	void <- write(c(bin_bash, cmd, cmd_keep_rds_only), bash_job_file)
}



#'  Write a batch script to merge rds files from filterVCF runs slurm job
#' 
#' This function should not be used directly but is part of the filterVCF.R
#' pipeline
#' 
#' @inheritParams write_batch_script_slurm
#' @return Write batch job files
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~batchCreator
#' @examples
#' 
#' 
#' @export
write_batch_script_slurm_merged_run_rds <- function(vcf_file, outputdir,...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}
	bin_bash <- sprintf("#!/bin/bash -l
#SBATCH -A snic2017-7-175
#SBATCH -J filter_vcf
#SBATCH -p core -n %s
#SBATCH -t %s

cd %s
", opt$no_cores,  "60:00:00" , dirname(outputdir))
	#cmd <- "Rscript test.R --missing 0.2 -I a --gzip_program gunzip"
	cmd <- get_call()
    # cmd <- call
	cmd <- gsub("--vcf_file .+vcf??? -|-I .+ -", paste0("--vcf_file " , basename(vcf_file), " -"), cmd)
	cmd_keep_rds_only <- ""

	cmd <- gsub("--dry_run", "", cmd)
	cmd <- paste(cmd, " --merge_only")
	message("setting --no_cores to 1, multithreading is not yet optimal, and sometimes create memmory allocation problem. Change it in the batch script to try multicore")
	cmd <- gsub("--no_cores .+ -", "--no_cores '1' -", cmd)

	bash_job_file <- sprintf("merge_run_rds_%s.sh", gsub(".vcf.gz|.vcf", "", basename(vcf_file)))
	void <- write(c(bin_bash, cmd, cmd_keep_rds_only), bash_job_file)
}



#'  Write a batch script to merge vcf and bed files from filterVCF runs slurm job
#' 
#' This function should not be used directly but is part of the filterVCF.R
#' pipeline
#' 
#' @inheritParams write_batch_script_slurm
#' @return Write batch job files
#' @author  ~~Benjamin Laenen~~
#' @seealso  objects to See Also as \code{\link{help}}, 
#' @references  
#' @keywords ~batchCreator
#' @examples
#'
#' 
#' @export
write_batch_script_slurm_merged_run_chr <- function(split_VCF_filename, outputdir, bed_output_dir, ...){
	myDots <- list(...)
	if (!is.null(myDots$opt)){
		opt <- myDots$opt
	}else{
		opt <- parse_args(OptionParser(option_list=initialise_option()))
	}
	bin_bash <- sprintf("#!/bin/bash -l
#SBATCH -A snic2017-7-175
#SBATCH -J merged_run_chr
#SBATCH -p core -n %s
#SBATCH -t %s

cd %s
", opt$no_cores,  "60:00:00" , escape_dropbox_pro_path(outputdir))

	cmd <- ""

	#merging bed files
	outputdirs <- gsub(".vcf.gz", "", split_VCF_filename)
	bed_files <- lapply(outputdirs, list.files, pattern = ".bed", recursive = TRUE, full.names = TRUE)

	#we can find a more clever way by checking the names
	#of the bedfile that meatch each pther except the chr
	#for now we considered that they are equal and merge
	#each by position in the list
	bed_files <- do.call("cbind", bed_files)

	filter_name <- gsub(gsub(".vcf.gz|.vcf", "", basename(opt$vcf_file)), "", basename(bed_files[,1]))
	filter_name <- sapply(sapply(strsplit(filter_name, "_"), "[", -c(1:2)), paste, collapse = "_")

	bed_files <- cbind(filter_name, bed_files)

	#keeping the import bed file in the output directory
	cmd_merge_bed <- apply(bed_files, 1, function(x) sprintf("cat %s > %s/%s", paste(sapply(x[-1], escape_dropbox_pro_path), collapse = " "), escape_dropbox_pro_path(bed_outputdir),paste0(gsub(".vcf.gz|.vcf", "", basename(opt$vcf_file)), "_", x[1])))
	cmd_merge_bed <- gsub("(^.+/)bed_filters/(.+_invariable.bed)$", "\\1\\2", cmd_merge_bed)
	cmd_merge_bed <- gsub("(^.+/)bed_filters/(.+_removed_sites.bed)$", "\\1\\2", cmd_merge_bed)

	#merging the vcf file
	vcf_files <- lapply(outputdirs, list.files, pattern = ".vcf.gz", recursive = TRUE, full.names = TRUE)


	cmd_header <- sprintf("zgrep '^#' %s > header", escape_dropbox_pro_path(vcf_files[1]))
	cmd_vcf <- sprintf("zgrep -v '^#' %s > %s.vcf.gz", paste(sapply(vcf_files, escape_dropbox_pro_path), collapse=" "), paste0(gsub(".vcf.gz|.vcf", "", basename(opt$vcf_file)), "_", "filtered"))

	#removing the bedfiles and vcffiles by chr?

	bash_job_file <- sprintf("merge_run_by_chr_%s.sh", gsub(".vcf.gz|.vcf", "", basename(opt$vcf_file)))
	void <- write(c(bin_bash, cmd_merge_bed, cmd_header, cmd_vcf), bash_job_file)

}
