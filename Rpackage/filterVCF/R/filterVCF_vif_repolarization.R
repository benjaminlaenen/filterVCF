#' Create input to run P.Keightley 'est-sfs'
#' 
#' This function takes a vif object and create an input to run P.Keightly
#' program that aims at estimating the unfolded SFS using two outgroups in a
#' maximum likelihood framework.
#' It is expected that the slot @divergence of the vif object is 
#' filled and two outgroups are present.
#' The program est-sfs require that all sites have the same number of
#' alleles. If we don't want to loose too much, we can approximate the
#' frequency with freq / (nb sampled allele * max nb of sampled alleles). This
#' approxiation is robust if there is not too much missing data.
#' 
#' Below is the descriptoin from the manual on the input specifications:
#'	# The format is the following:
#'	#     There are 12 space-separated columns, 4 for the focal species and 4
#'	# for each of the two outgroups. If the single outgroup program is used,
#'	# the last 4 columns are ignored, but must be present (they can just be
#'	# a duplication of outgroup 1). Each line in the input file gives the
#'	# count of alleles at a site. The order of bases is A C G T.
#'
#'	# For example, the first line in the example data file is:
#'
#'	# 20 0 0 0                0 0 0 1 0 0 0 1
#'
#'	# At this site, all of the n = 20 copies sampled are A in the focal
#'	# species. In the outgroups, a single copy has been sampled, and both
#'	# are T. All sites must have the same number of copies sampled in the
#'	# focal species and one copy sampled in each of the outgroups.
#' 
#' 
#' @param vif a VIF object
#' @param approximate_freq Should the frequency by approximated as explained above?
#' @return A GRange object with metadata suitable to create an input for est-sfs. The values associated with each sitesare the following
#' \item{freq }{allele frequency of the alternates alleles} %% 
#' \item{nsample }{total number of chromosomes sampled} %% ...
#' \item{REF }{base of the REF} %%
#' \item{ALT }{base of the ALT} %%
#' \item{OUT1 }{base of the first outgroup} %%
#' \item{OUT2 }{base of the second outgroup} %%
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~est-sfs
#' @examples
#' 
#' input_est-sfs <- make_GRange_freq_base_Keightley(vif_object)
#' input_est-sfs
#' 
#' @export
make_GRange_freq_base_Keightley <- function(vif, approximate_freq = TRUE,...){
	# The format is the following:
	#     There are 12 space-separated columns, 4 for the focal species and 4
	# for each of the two outgroups. If the single outgroup program is used,
	# the last 4 columns are ignored, but must be present (they can just be
	# a duplication of outgroup 1). Each line in the input file gives the
	# count of alleles at a site. The order of bases is A C G T.

	# For example, the first line in the example data file is:

	# 20 0 0 0                0 0 0 1 0 0 0 1

	# At this site, all of the n = 20 copies sampled are A in the focal
	# species. In the outgroups, a single copy has been sampled, and both
	# are T. All sites must have the same number of copies sampled in the
	# focal species and one copy sampled in each of the outgroups.
	get_table_ACGT <- function(all_site_freq_in_div){
	   nsample <- max(elementMetadata(all_site_freq_in_div)["nsample"][,1])
	   base = c("A", "C", "G", "T")
	
	   ACGT <- matrix(0, nrow = length(all_site_freq_in_div), ncol = 4)
	   # adding invariant
	   index_invariant <- elementMetadata(all_site_freq_in_div)["ALT"][,1] == "*"
	   for(i in seq_along(base)){
		ACGT[index_invariant &    elementMetadata(all_site_freq_in_div)["REF"][,1] == base[i] , i] <- nsample
	   }
	
	   # add snp and fixed
	   index_snp <- elementMetadata(all_site_freq_in_div)["ALT"][,1] != "*"
	   for(i in seq_along(base)){
		ACGT[index_snp &    elementMetadata(all_site_freq_in_div)["REF"][,1] == base[i] , i] <- nsample - elementMetadata(all_site_freq_in_div)[index_snp & elementMetadata(all_site_freq_in_div)["REF"][,1] == base[i], "freq"]
		ACGT[index_snp &    elementMetadata(all_site_freq_in_div)["ALT"][,1] == base[i] , i] <-           elementMetadata(all_site_freq_in_div)[index_snp & elementMetadata(all_site_freq_in_div)["ALT"][,1] == base[i], "freq"]
	   }
	   elementMetadata(all_site_freq_in_div) <- ACGT
	   return(all_site_freq_in_div)
	}

	get_table_ACGT_div <- function(all_site_div_in_freq){
		nsample <- 1
		#check case of the divergence
		#only checking a sample
		bases_div <- unique(unlist(sapply(1:ncol(elementMetadata(all_site_div_in_freq)), function(i) as.character(unique(elementMetadata(all_site_div_in_freq)[,i])))))
		case_div <- check_case(bases_div, isupper=TRUE)
		if(all(case_div) | all(!case_div)){
			if(all(case_div)){
				base = c("A", "C", "G", "T")
			}else{
				base = c("a", "c", "g", "t")
			}
		}else{
			message("Sites in the divergence have mixed cases, check the input divergence. Changing all to upper case.")
			for(i in seq_len(ncol(elementMetadata(all_site_div_in_freq)))){
				elementMetadata(all_site_div_in_freq)[,i] <- toupper(elementMetadata(all_site_div_in_freq)[,i])
			}
		}
		ACGT_div <- matrix(0, nrow = length(all_site_div_in_freq), ncol = 8)
	
		for(i in seq_along(base)){
			ACGT_div[elementMetadata(all_site_div_in_freq)[,2] == base[i] , i] <- 1
		}
	
		if(ncol(elementMetadata(all_site_div_in_freq)) > 2){
			for(i in seq_along(base)){
				ACGT_div[elementMetadata(all_site_div_in_freq)[,3] == base[i] , i+4] <- 1
			}
		}else{
			ACGT_div[, 5:8] <- ACGT_div[, 1:4]
		}

		elementMetadata(all_site_div_in_freq) <- ACGT_div
		return(all_site_div_in_freq)
	}

	# Create a Granges with frequency
	# and nb of missing for all sites
	# Note that it is assumed that the fixed
	# and invariant are complete.
	# In other word, we assume that if
	# we were to sample the missing genotype
	# in the invariant, they will remain invariant.
	# This assumption holds if the sample used to discover the snp
	# is large enough and missing threshold not too low.    
	dots <- list(...)
	if(!is.null(dots$use_ff)){
		use_ff <- dots$use_ff
	}else{
		use_ff <- FALSE
	}
	freq <- get_frequency(vif, use_ff = use_ff)
	nsample <- max(elementMetadata(freq$vcf_freq)["nsample"][,1])
	
	#check that fixed sites have a width of 1
	#if not it means that there is some interval
	#and the replacement wont work to assign 
	#ALT base and outgroup on a per site basis.
	if(any(unique(width(freq$fixed)) != 1)){
		freq$fixed <- expand_all_pos_Granges(freq$fixed)
	}

	elementMetadata(freq$invariant)["freq"] <- 0
	elementMetadata(freq$fixed)["freq"] <- nsample
	elementMetadata(freq$invariant)["nsample"] <- elementMetadata(freq$fixed)["nsample"] <- nsample
	# Getting the alternate alleles
	# to produce the input file with ACGT columns
	elementMetadata(freq$invariant)["ALT"] <- "*"
	elementMetadata(freq$vcf_freq)["ALT"]  <- elementMetadata(intersect_Granges(alternate(vif), freq$vcf_freq, merge_interval=FALSE))["ALT"]
	elementMetadata(freq$fixed)["ALT"]     <- elementMetadata(intersect_Granges(alternate(vif), freq$fixed, merge_interval=FALSE))["ALT"]

	all_site_freq <- sort(c(freq$vcf_freq, freq$invariant, freq$fixed))

	# keep only sites present in the divergence
	# as we will only be able to determinae the
	# ancestral state for those sites.
	all_site_freq_in_div <- intersect_Granges(all_site_freq, divergence(vif), merge_interval=FALSE)

	#get the base from the reference genome for each positions
	elementMetadata(all_site_freq_in_div)["REF"] <- as.character(getSeq(reference(vif), all_site_freq_in_div))

	if(approximate_freq){
		# The program est-sfs require that all a 
		# sites have the same number of alleles
		# If we don't want to loose too much 
		# I approximate the freq but upsampling
		# as freq / n smapled allele * max nb of sampled alleles.
		# This approxiation is robust if there is not
		# too much missing data.
		elementMetadata(all_site_freq_in_div)["freq"] <- floor(elementMetadata(all_site_freq_in_div)["freq"][,1] / elementMetadata(all_site_freq_in_div)["nsample"][,1] * nsample)
	}else{
		# Remove sites with missing
		index_missing <- elementMetadata(all_site_freq_in_div)["nsample"][,1] != nsample
		all_site_freq_in_div <- all_site_freq_in_div[!index_missing]
	}

	ACGT_freq <- get_table_ACGT(all_site_freq_in_div)

	# get the outgroup part
	# ACGT ACGT
	all_site_div_in_freq <- intersect_Granges(divergence(vif), all_site_freq, merge_interval=FALSE)
	ACGT_div <- get_table_ACGT_div(all_site_div_in_freq)

	ACGT <- ACGT_freq
	# add the column for the outgroups
	elementMetadata(ACGT)[,5:12] <- elementMetadata(ACGT_div)

	#test if uppercase and change to upper if not
	# use 10000 random sample to test
	# it wont check all sites
	elementMetadata(ACGT)["REF"] <- elementMetadata(all_site_freq_in_div)["REF"]
	elementMetadata(ACGT)["ALT"] <- elementMetadata(all_site_freq_in_div)["ALT"]
	elementMetadata(ACGT)["OUT1"] <- elementMetadata(all_site_div_in_freq)[,2]
	elementMetadata(ACGT)["OUT2"] <- elementMetadata(all_site_div_in_freq)[,3]

	if(any(check_case(elementMetadata(ACGT)["REF"][,1], Sample = 10000, isupper=FALSE))){
		elementMetadata(ACGT)["REF"] <- toupper(elementMetadata(ACGT)["REF"][,1])
	}
	if(any(check_case(elementMetadata(ACGT)["ALT"][,1], Sample = 10000, isupper=FALSE))){
		elementMetadata(ACGT)["ALT"] <- toupper(elementMetadata(ACGT)["ALT"][,1])
	}
	if(any(check_case(elementMetadata(ACGT)["OUT1"][,1], Sample = 10000, isupper=FALSE))){
		elementMetadata(ACGT)["OUT1"] <- toupper(elementMetadata(ACGT)["OUT1"][,1])
	}
	if(any(check_case(elementMetadata(ACGT)["OUT2"][,1], Sample = 10000, isupper=FALSE))){
		elementMetadata(ACGT)["OUT2"] <- toupper(elementMetadata(ACGT)["OUT2"][,1])
	}
	return(ACGT)
}


#' Write input for est-sfs in a file or internally
#' 
#' This function create a caracter vector using the output from
#' make_GRange_freq_base_Keightley() to write it to a file if an output file is
#' specified or return a character with the input.
#' 
#' 
#' @param ACGT GRanges output from make_GRange_freq_base_Keightley()
#' @param output output file name, if NULL the function only return a character vctor with the input to run est-sfs
#' @return a character vector with the input to run est-sfs
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~est-sfs
#' @examples
#' 
#' input_est-sfs <- make_GRange_freq_base_Keightley(vif_object)
#' write_input_Keightley(input_est-sfs)
#' 
#' 
#' @export
write_input_Keightley <- function(ACGT, output=NULL, ...){
	# output = "out_Keightley"
	#tmp <- tempfile()

	input <- as.data.frame(elementMetadata(ACGT)[,1:12])

	input_Keightley <- paste(paste(input[,1], input[,2], input[,3], input[,4], sep = " "),
						     paste(input[,5], input[,6], input[,7], input[,8], sep = " "),
						     paste(input[,9], input[,10], input[,11], input[,12], sep = " "), sep = "\t")
	if(!is.null(output)){
		if(file.exists(output)) file.remove(output)
			writeLines(input_Keightley, output)
		}
	return(input_Keightley)
}


#' Wrapper to run est-sfs from P.Keightley 
#' 
#' This function runs est-sfs from P.Keightley. The program need to be
#' installed and and be callable using a system() call. This has only been test
#' on UNIX (mac and linux). This is used internally in the larger function
#' make_repol_key() that should be used to repolarize an vif or filterVCF
#' object. However, this wrapper function allows more freedom in the model.
#' 
#' The aim of the function is to return the file with the ancestral
#' probability for each sites but full output from the program can also be
#' available if write_files = TRUE
#' 
#' 
#' @param ACGT output from make_GRange_freq_base_Keightley()
#' @param Model = "JC", "kimura", "rate6". Sequence evolution model to use. A LRT can be conducted ad hoc to shoose the best model
#' @param nrandom = 5 number of random start
#' @param noutgroups = 2 Number of outgroup (fixed at 2 for the moment)
#' @param bin_path = "est-sfs" location of the executable
#' @param write_files=FALSE Should a full output be written to a folder
#' @param name = "" common root name for output
#' @param ... options (NOT USED)
#' @return A data.frame with the same number of sites as the input Grange. It contains the probability that the major allele is ancestral (see est-sfs manual for more detail)
#' @author  ~~Benjamin Laenen~~
#' @keywords ~est-sfs
#' @examples
#' 
#' input_est-sfs <- make_GRange_freq_base_Keightley(vif_object)
#' run_Keightley_estsfs(input_est-sfs, Model = "kimura", bin_path = "est-sfs", write_files = TRUE)
#' 
#' @export
run_Keightley_estsfs <- function(ACGT, Model = c("JC", "kimura", "rate6"), nrandom = 5, noutgroups = 2, bin_path = "est-sfs", write_files=FALSE, name = "", ...){
    input_Keightley <- write_input_Keightley(ACGT)
    if(length(Model) == 3) Model = "JC"
    if(Model == "JC") model = 0
    if(Model == "kimura") model = 1
    if(Model == "rate6") model = 2
    
    tmpdir <- paste("ESTSFS", round(runif(1,1000,90000)), sep = "_")
    dir.create(tmpdir)
    # start_dir    <- getwd()
    
    # setwd(tmpdir)
    #create the config file
    config <- rbind(n_outgroup = noutgroups, model = model, nrandom = nrandom)
    config_file <- paste(tmpdir, paste("config", Model, sep ="_"), sep = .Platform$file.sep)
    write.table(config, config_file, col.names = FALSE, quote = FALSE)

    #create the input
    input_file <- paste(tmpdir, paste("input_Keightley", name, sep ="_"), sep = .Platform$file.sep)
    writeLines(input_Keightley, input_file)
 
    #create the seed file
    seed_file <- paste(tmpdir, "seed", sep = .Platform$file.sep)
    writeLines(as.character(round(runif(1, 1000, 10000000))), seed_file)

    output_sfs <- paste(tmpdir,paste("ouput", Model, "sfs", name, sep ="_"), sep = .Platform$file.sep)
    output_p_anc <- paste(tmpdir,paste("ouput", Model, "p_anc", name, sep ="_"), sep = .Platform$file.sep)

    #create the command to run
    cmd <- paste(escape_dropbox_pro_path(bin_path), config_file, input_file, seed_file, output_sfs, output_p_anc)
    #run the program
    system(cmd)
    
    # find the header (varying number of line
    # depending on the model)
    skip <- readLines(output_p_anc, n =20)
    skip <- max(grep("^0", skip))
    p_anc <- read.table(output_p_anc, skip = skip)[,-2]
    colnames(p_anc) <- c("sites_nb", "p_major_ancestral", "p_node1_A", "p_node1_C", "p_node1_G", "p_node1_T")

    # setwd(start_dir)
    if(write_files){
		void <- file.copy(paste(tmpdir, basename(output_p_anc), sep = "/"), ".", overwrite = TRUE)
		void <- file.copy(paste(tmpdir, basename(input_file), sep = "/"), ".", overwrite = TRUE)
		void <- file.copy(paste(tmpdir, basename(output_sfs), sep = "/"), ".", overwrite = TRUE)
    }
    unlink(tmpdir, recursive = TRUE)
    return(p_anc)
}

#' Transform p_anc to repolKey
#' 
#' Internal function used in make_repol_key(). Only experienced users should use this function.
#' 
#' 
#' @param p_anc output from run_Keightley_estsfs()
#' @param ACGT output from make_GRange_freq_base_Keightley()
#' @param pp_threshold The cutoff over which the probability of ancestral is consider. For example, if pp_threshold = 0.8, sites with probability  <0.2 and >0.8 will be considered for repolarization. 
#' @return a repolKey with information to repolarize the vif.
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~est-sfs
#' @examples
#' 
#' input_est-sfs <- make_GRange_freq_base_Keightley(vif_object)
#' p_anc <- run_Keightley_estsfs(input_est-sfs, Model = "kimura", bin_path = "est-sfs", write_files = TRUE)
#' repolKey(vif_object) <- p_anc2repol_key(p_anc, input_est-sfs , pp_threshold = 0.8)
#' 
#' @export
p_anc2repol_key <- function(p_anc, ACGT, pp_threshold = 0.8, ... ){
    # Function    that read the p_anc file from est-sfs
    # and create a Granges object with the probability that
    # that the REF is ancestral
    
	# Test if p_anc and ACGT have the same length
	if(nrow(p_anc) != length(ACGT)){
		message(sprintf("Output length from est-sfs differ from input : %s / %s", nrow(p_anc) , length(ACGT)))
		quit("no")
	}
    elementMetadata(ACGT)[,"p_major_ancestral"] <- p_anc$p_major_ancestral
    elementMetadata(ACGT)[,"p_node1_A"] <- p_anc$p_node1_A
    elementMetadata(ACGT)[,"p_node1_C"] <- p_anc$p_node1_C
    elementMetadata(ACGT)[,"p_node1_G"] <- p_anc$p_node1_G
    elementMetadata(ACGT)[,"p_node1_T"] <- p_anc$p_node1_T
    #there is an issue with the probability of
    #major being ancestral for invariant.
    #The p*_anc is always 1 even if the outgroup
    #is different. I will use the probabilty of the
    #ancestral base at the node1 ((in,ou1)node1, out2)
    #tp infer which base is ancestral and if the sites must
    #be repolarized.
    #Several case are to be considered
    # 1) in != out1 != out2 : high uncertainty in ancestral. either use the highest prob or remove the sites
    # 2) 
    # 
    index_max <- apply(as.data.frame(elementMetadata(ACGT)[,c("p_node1_A", "p_node1_C", "p_node1_G", "p_node1_T")]),1, which.max)
    anc_base <- index_max
    anc_base[anc_base==1] <- "A"
    anc_base[anc_base==2] <- "C"
    anc_base[anc_base==3] <- "G"
    anc_base[anc_base==4] <- "T"
    elementMetadata(ACGT)[,"anc_base"] <- anc_base
    elementMetadata(ACGT)[index_max==1,"anc_base_prob"] <- elementMetadata(ACGT)[index_max==1,"p_node1_A"]
    elementMetadata(ACGT)[index_max==2,"anc_base_prob"] <- elementMetadata(ACGT)[index_max==2,"p_node1_C"]
    elementMetadata(ACGT)[index_max==3,"anc_base_prob"] <- elementMetadata(ACGT)[index_max==3,"p_node1_G"]
    elementMetadata(ACGT)[index_max==4,"anc_base_prob"] <- elementMetadata(ACGT)[index_max==4,"p_node1_T"]

    index_major <-apply(as.data.frame(elementMetadata(ACGT)[,1:4]),1, which.max)
    major_base <- index_major
    major_base[major_base==1] <- "A"
    major_base[major_base==2] <- "C"
    major_base[major_base==3] <- "G"
    major_base[major_base==4] <- "T"
    elementMetadata(ACGT)[,"major"] <- major_base

    #index inv
    index_inv <- elementMetadata(ACGT)[,"ALT"] == "*"

    #test if ancestral is different from REF
    index_ref_not_anc <- elementMetadata(ACGT)[,"anc_base"] != elementMetadata(ACGT)[,"REF"]

    # elementMetadata(ACGT)[index_ref,] # do not repol
    # elementMetadata(ACGT)[index_alt,] # repol according to p_major_ancestral > 0.8
    # elementMetadata(ACGT)[index_notrefnotalt,] # repol according to p_major_ancestral > 0.8 with major as ancestral if >0.8 and minor if p<0.2
    # elementMetadata(ACGT)[index_notref_inv,] #repol according to anc_base

    elementMetadata(ACGT)[,"repolarize_ml"] <- 0

    index_major_alt_is_anc <- elementMetadata(ACGT)[,"ALT"] == elementMetadata(ACGT)[,"major"] & elementMetadata(ACGT)[,"p_major_ancestral"] > pp_threshold
    index_major_ref_not_anc <- elementMetadata(ACGT)[,"REF"] == elementMetadata(ACGT)[,"major"] & elementMetadata(ACGT)[,"p_major_ancestral"] < 1-pp_threshold
    
    elementMetadata(ACGT)[!index_inv & index_major_alt_is_anc,"repolarize_ml"] <- 1
    elementMetadata(ACGT)[!index_inv & index_major_ref_not_anc,"repolarize_ml"] <- 1
    # correct for the fact that est-sfs
    # assign a prob of major anc = 1
    # to invariant even when the ancestral base is
    # different from the ref.
    # Here I assume that the ancestral base
    # is the one recovered from the model.
    elementMetadata(ACGT)[index_inv & index_ref_not_anc,"repolarize_ml"] <- 1

    #site that we cannot call
    #will be missing from the key
    #but we could use the est-sfs 
    #later without using the key
    #to compare the results of the usfs
    index_uncertain_prob_major_anc <- elementMetadata(ACGT)[,"p_major_ancestral"] < pp_threshold & elementMetadata(ACGT)[,"p_major_ancestral"] > 1-pp_threshold
    elementMetadata(ACGT)[index_uncertain_prob_major_anc,"repolarize_ml"] <- 2

    #case of the 0.5 frequency
    #which allele is major?
    #do not repol if anc_base == REF
    #repol if anc_base == ALT
    #NA if anc_base != ALt or REF
    index_max_0.5freq <-apply(as.data.frame(elementMetadata(ACGT)[,1:4]),1, function(x) max(x) == sum(x)/2)
    index_alt <- elementMetadata(ACGT)[,"anc_base"] == elementMetadata(ACGT)[,"ALT"]
    index_ref <- elementMetadata(ACGT)[,"anc_base"] == elementMetadata(ACGT)[,"REF"] 
    elementMetadata(ACGT)[index_max_0.5freq & index_ref,"repolarize_ml"] <- 0
    elementMetadata(ACGT)[index_max_0.5freq & index_alt,"repolarize_ml"] <- 1
    elementMetadata(ACGT)[index_max_0.5freq & !index_ref & !index_alt,"repolarize_ml"] <- 2

    res <- ACGT
    elementMetadata(res) <- elementMetadata(ACGT)[,c("repolarize_ml", "anc_base", "anc_base_prob", "major", "p_major_ancestral", "REF", "ALT", "OUT1", "OUT2")]
    return(res)
}

#' Create a repolKey using parsimony
#' 
#' Internal function used in make_repol_key(). Only experienced users should use this function.
#' 
#' 
#' @param ACGT output from make_GRange_freq_base_Keightley()
#' @param out1_out2_equal Condition that the two outgroups shouldhave the same base.
#' @return a repolKey with information to repolarize the vif.
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~est-sfs
#' @examples
#' 
#' input_est-sfs <- make_GRange_freq_base_Keightley(vif_object)
#' repolKey(vif_object) <- parsimony_repol_key(input_est-sfs , out1_out2_equal = TRUE)
#' 
#' @export
parsimony_repol_key <- function(ACGT, out1_out2_equal = TRUE){
	repolKey <- ACGT
	elementMetadata(repolKey)[,"repolarize_pars"] <- 0

    index_ref_out <- elementMetadata(ACGT)[,"REF"] == elementMetadata(ACGT)[,"OUT1"]
    index_alt_out <- elementMetadata(ACGT)[,"ALT"] == elementMetadata(ACGT)[,"OUT1"]
    index_alt_inv <- elementMetadata(ACGT)[,"ALT"] == "*"
    index_out1_out2_same <- elementMetadata(ACGT)[,"OUT2"] == elementMetadata(ACGT)[,"OUT1"]
    if(out1_out2_equal){
	    elementMetadata(repolKey)[index_out1_out2_same & !index_ref_out,"repolarize_pars"] <- 1
	    elementMetadata(repolKey)[index_out1_out2_same & !index_ref_out & !index_alt_out & index_alt_inv,"repolarize_pars"] <- 1
	    elementMetadata(repolKey)[!index_out1_out2_same,"repolarize_pars"] <- 2
    }else{
    	#using only the first outgroup
    	#that need to either match the ALT
    	#or different from the REF if ALT is *
    	#missing is ALT and REF are different from OUT1
	    elementMetadata(repolKey)[index_alt_out,"repolarize_pars"] <- 1
	    elementMetadata(repolKey)[!index_ref_out & !index_alt_out & index_alt_inv,"repolarize_pars"] <- 1
	    elementMetadata(repolKey)[!index_ref_out & !index_alt_out & !index_alt_inv,"repolarize_pars"] <- 2
    }

    elementMetadata(repolKey) <- elementMetadata(repolKey)[,c("repolarize_pars", "REF", "ALT", "OUT1", "OUT2")]
    return(repolKey)

    # utility to check if the key is working as expected
	# index_repol <- elementMetadata(repolKey)[,"repolarize"] == 1
    # as.data.frame(elementMetadata(repolKey)[index_repol,])[1:100,]
}


#' Create a repolarization key using parsimony or ML
#' 
#' This function uses divergence to estimates the ancestral states for each
#' sites that has been sequenced and is present in the vif object. The maximum
#' likelihood uses the external program est-sfs (see ?run_Keightley_estsfs() for
#' more detail). This function DOES NOT repolarize but just creates the repolKey.
#' 
#' 
#' @param vif a VIF object
#' @param methods "ML" "parsimony". Which methods to use to create the key? If both the repolKey will contains a column "repolarize_pars" and "repolarize_ML" which indicates which sitesshould be switch when repolarizing.
#' @return a repolKey
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~repolarization
#' @examples
#' 
#' repolKey(vif_object) <- make_repol_key(vif_object, methods = c("ML", "parsimony"))
#' 
#' 
#' @export
make_repol_key <- function(vif, methods = c("ML", "parsimony"), ...){
	ACGT <- make_GRange_freq_base_Keightley(vif, approximate_freq=TRUE, ...)
	if(any(methods %in% "ML")){
		p_anc <- tryCatch(run_Keightley_estsfs(ACGT, Model = "kimura"), error=function(e) NULL)
		# if the file for input in est-sfs is too large (did not work for file
		# > 1Gb) the program cannot run because of memmory allocation and
		# return a core dumped. We can split the ACGT input into chromomose
		# and run the est-sfs- program again
		#139 mean segmentation faults
		if(is.null(p_anc)){
			chr <- levels(ACGT@seqnames)
			
			# if the number of core is not specified
			# take the number of chr as no_cores
			if(is.null(list(...)$opt$no_cores)){
				core_present <- detectCores()
				if(length(chr) > core_present){
					opt <- list(no_cores = core_present)
				}else{
					opt <- list(no_cores = length(chr))
				}
			}

			#split by chromosome and run on multicore
			message("input is too large for est-sfs, splitting by chromosome and running using multicore.")
			ACGT_by_chr <- split(ACGT, ACGT@seqnames)
			p_anc <- mclapply_try(seq_along(chr), function(i) run_Keightley_estsfs(ACGT_by_chr[[i]], Model = "kimura", name = chr[i]), mc.cores = opt$no_cores)
			
			# if mclapply fails, try to run lapply on element that are NULL
			# Sometimes multicore fails for some reason but only NULL
			# is reported without warnings.
			if(any(sapply(p_anc, is.null))){
				# removing the tmp folder if the run failed
				failed_run <- list.files(pattern = "ESTSFS_", include.dirs = TRUE, full.names = TRUE)
				void <- unlink(failed_run, recursive = TRUE)
				p_anc_NULL <- lapply(seq_along(chr)[sapply(p_anc, is.null)], function(i) run_Keightley_estsfs(ACGT_by_chr[[i]], Model = "kimura", name = chr[i]))
				p_anc[sapply(p_anc, is.null)] <- p_anc_NULL
			}
			p_anc <- do.call("rbind", p_anc)

		}

		repolKey_ML <- p_anc2repol_key(p_anc, ACGT)
		repolKey <- repolKey_ML
	}
	if(any(methods %in% "parsimony")){
		repolKey_parsimony <- parsimony_repol_key(ACGT, out1_out2_equal = TRUE)
		repolKey <- repolKey_parsimony
	}
	if (all(methods == c("ML", "parsimony"))){
		repolKey <- repolKey_ML
		elementMetadata(repolKey)[,"repolarize_pars"] <- elementMetadata(repolKey_parsimony)[,"repolarize_pars"]
	}
	return(repolKey)

	# way to check how much parsimony and ml differs
	# table(as.data.frame(elementMetadata(repolKey)[,c("repolarize_ml", "repolarize_pars")]))
}


#' Repolarize a vif_object
#' 
#' Use the information in the repolKey to polarize the SNPs, the invariant and
#' the fixed in the vif object. Note that the repolKey must be created
#' beforehand.
#' 
#' @param vif a VIF object
#' @param methods "ML" OR "parsimony". Which methods to use to repolarize? Note that the repolKey must have been create with the methods in ordor to use it (logical!)
#' @return a polarized vif object
#' 
#' @author  ~~Benjamin Laenen~~
#' @keywords ~repolarization
#' @examples
#' 
#' repolKey(vif_object) <- make_repol_key(vif_object, methods = c("ML", "parsimony"))
#' polarized_vif  <- repolarize_VIF(vif_object, methods = "ML")
#' 
#' 
#' @export
repolarize_VIF <- function(vif, methods = "ML"){
	index_repol <- elementMetadata(repolKey(vif))[,"repolarize_ml"] == 1
	index_drop_repol <- elementMetadata(repolKey(vif))[,"repolarize_ml"] == 2
	
	res <- VIF()
	#start by repolarizing invariant and 
	#drop site that we cannot asses by ML
	repol_fixed <- intersect_Granges(invariant(vif), repolKey(vif)[index_repol])  
	repol_inv <- intersect_Granges(fixed(vif), repolKey(vif)[index_repol])  
	
	invariant(res) <- c(repol_inv, intersect_Granges(invariant(vif), repol_fixed, invert = TRUE))
	fixed(res)     <- c(repol_fixed, intersect_Granges(fixed(vif), repol_inv, invert = TRUE))

	invariant(res) <- intersect_Granges(invariant(res), repolKey(vif)[index_drop_repol], invert = TRUE)
	fixed(res)     <- intersect_Granges(fixed(res), repolKey(vif)[index_drop_repol], invert = TRUE)
	
	#change the alternates states
	alternate(res) <- intersect_Granges(alternate(vif), repolKey(vif)[index_repol], merge_interval = FALSE)
	#change the column names to shift the REF to ALT
	#and vice versa then reorder and add the sites non repolarized.
	colnames(elementMetadata(alternate(res))) <- c("ALT", "REF")
	elementMetadata(alternate(res)) <- elementMetadata(alternate(res))[,c("REF", "ALT")]

	alternate(res) <- c(alternate(res), intersect_Granges(alternate(vif), repolKey(vif)[!index_repol & !index_drop_repol], merge_interval = FALSE))
	alternate(res) <- sort(alternate(res))
	
	#change the VCF state 0/0 to 1/1
	#and the REF aLT columns
	vcf(res) <- repolarize_vcf(vcf(vif), repolKey(vif), methods = "ML")

	#add the sites that cannot be plarized
	# to the removed sites
	removed_sites(res) <- c(removed_sites(res), vcf2Grange(intersect_Granges(vcf(vif), repolKey(vif)[index_drop_repol])))
	removed_sites_inv(res) <- c(removed_sites_inv(res), intersect_Granges(vcf(vif), repolKey(vif)[index_drop_repol]))
	removed_sites_inv(res) <- c(removed_sites_inv(res), intersect_Granges(fixed(vif), repolKey(vif)[index_drop_repol]))

	GetOpt(res) <- GetOpt(vif)	
	repolKey(res) <- repolKey(vif)
	divergence(res) <- divergence(vif)
	reference(res) <- reference(vif)
	res@bootstrap <- vif@bootstrap
	return(res)
}

#' Polarize a vcfR
#' 
#' Internal function used in repolarize_VIF(). Only experienced users. NOte that invariant and fixed will not be polarized here.
#' 
#' 
#' @param vcf vcfR object
#' @param repolKey a repolarization key create using make_repol_key()
#' @param methods "ML" OR "parsimony". Which methods to use to repolarize? Note that the repolKey must have been create with the methods in ordor to use it (logical!)
#' @return a polarized vcfR object
#' 
#' @author  ~~Benjamin Laenen~~
#' @references  put references to the literature/web site here
#' @keywords ~repolarization
#' @examples
#' 
#' repolKey(vif_object) <- make_repol_key(vif_object, methods = c("ML", "parsimony"))
#' polarized_VCF  <- repolarize_VIF(vcf(vif_object), methods = "ML")
#' 
#' 
#' @export
repolarize_vcf <- function(vcf, repolKey, methods = "ML"){
	if(methods == "ML") repol_column <- "repolarize_ml"
	if(methods == "parsimony") repol_column <- "repolarize_pars"
	index_drop_repol <- elementMetadata(repolKey)[,"repolarize_ml"] == 2
	original_length <- nrow(vcf)
	vcf <- intersect_Granges(vcf, repolKey[!index_drop_repol])
	message(sprintf("%s sites dropped with uncertainty in the repolarization", original_length - nrow(vcf)))
	repolKey_same_as_vcf <- intersect_Granges(repolKey[!index_drop_repol], vcf, merge_interval = FALSE)

	#change the ALT REF columna
	index_repol_vcf <- elementMetadata(repolKey_same_as_vcf)[,repol_column] == 1
	tmp_alt <- vcf@fix[index_repol_vcf,"ALT"] 
	vcf@fix[index_repol_vcf,"ALT"] <- vcf@fix[index_repol_vcf,"REF"]
	vcf@fix[index_repol_vcf,"REF"] <- tmp_alt
	rm(tmp_alt)

	# change the genotypes 0/0 -> 1/1
	# and 1/1 to 0/0 (or with | separator).
	# I dont change the hets 0/1 (quicker).
	# I use the package mgsub that can do
	# multiple search replacement without
	# intermediate objects.
	# Note that I dont change the other fields
	# AD or GQ, so we must use the info later!!
	
	vcf@gt[index_repol_vcf] <- mgsub(pattern = c("^0(.)0:", "^1(.)1:"), replacement = c("1\\11:", "0\\10:") , vcf@gt[index_repol_vcf])
	# adding the info that we changed the sites 
	# in the FILTER colmun
	vcf@fix[index_repol_vcf,"FILTER"] <- repol_column
	message(sprintf("%s sites repolarized", sum(index_repol_vcf)))
	return(vcf)
}

