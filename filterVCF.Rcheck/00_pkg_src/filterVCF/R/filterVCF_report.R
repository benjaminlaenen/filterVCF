#' Create a html report of filtration 
#' 
#' This function takes a filter VCF object that has been created by main_filter_VCF. 
#' 
#' The function is intended to be used during the pipeline filterVCF.R and
#' needs some slots to have a value, e.g. the filters and invariable.
#' 
#' @param RES a filterVCF object
#' @param outputdir Output directory to write the report
#' @param opt list of option used to filter the data. If NULL it is taken from the filterVCF@opt 
#' @return NULL  Write a report to a html with name taken from opt$report
#' @author  ~~Benjamin Laenen~~
#' @references  
#' @examples
#' 
#' 
#' RES <- main_filter_VCF(vcf_file, opt=opt)
#' create_report(RES,  outputdir = ".", opt = NULL)
#' 
#' 
#' @export
#' @import R2HTML
#' @import ape
create_report <- function(RES, outputdir = ".", opt = NULL, skip_chromqc = FALSE){
	if(is.null(opt)) opt <- GetOpt(RES)
	#Create a table of the number of sites each filter removed and how much data is variable, invariable, report also the sfs
	vcf_name <- basename(opt$vcf_file)
	#if(!exists("opt")) opt <- RES@opt
	make_plot_nb_sites <- function(RES, report_dir,  opt = NULL){
		if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))
		nb_sites <- list()
		nb_sites$total_sites <- nrow(vcfRaw(RES)) + sum(width(invariantRaw(RES))) + sum(width(fixed(RES)))
		nb_sites$empty <- ""
		nb_sites$nb_variable <- nrow(vcfRaw(RES))
		nb_sites$nb_variable_filtered <-  nrow(vcf(RES))
		nb_sites$nb_invariable <- sum(width(invariantRaw(RES)))
		nb_sites$nb_invariable_filtered <-  sum(width(invariant(RES)))
		nb_sites$plot_file <- paste0(report_dir, "/plots/", vcf_name, "_nb_sites_summary.pdf")

		pdf(nb_sites$plot_file)
		barplot(t(as.matrix(data.frame(c(nb_sites$total_sites,nb_sites$nb_variable,nb_sites$nb_invariable), c(0,nb_sites$nb_variable_filtered, nb_sites$nb_invariable_filtered)))), beside =TRUE, names.arg = c("Original", "Variable", "Invariable"), col = c("aquamarine2", 0, "darkolivegreen4", "darkolivegreen", "goldenrod3", "goldenrod1") )
		graphics.off()

		nb_sites$table_report <- as.data.frame(matrix(unlist(nb_sites[-7]), ncol =2, byrow =TRUE))
		colnames(nb_sites$table_report) <- c("unfiltered", "sites passing filters")
		rownames(nb_sites$table_report) <- c("Original", "Variable", "Invariable")
		nb_sites$plot_file <- paste0("./plots/", basename(nb_sites$plot_file))
		return(nb_sites)
	}

	create_plot_filter <- function(RES, report_dir, opt = NULL){
		if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))

		#nb of sites removed by each filter
		snp_remove_by_filter <- data.frame(snp_sites_removed=unlist(lapply(Filters(RES), sum, na.rm = TRUE)))
		snp_remove_by_filter$pc_total <- snp_remove_by_filter$snp_sites_removed / nrow(vcfRaw(RES))
		inv_remove_by_filter <- data.frame(inv_sites_removed=unlist(lapply(RES@filters_inv, sum, na.rm = TRUE)))
		inv_remove_by_filter$pc_total <- inv_remove_by_filter$inv_sites_removed / sum(width(invariantRaw(RES)))

		if(nrow(inv_remove_by_filter) == 0) inv_remove_by_filter <- data.frame(inv_sites_removed=c(0,0),pc_total=c(0,0))
		plot_file_snp <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_filters_variable.pdf")
		plot_file_inv <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_filter_invariable.pdf")

		pdf(plot_file_snp)
		par(mar = c(12, 3, 2, 2), cex.axis = 0.6, mfrow = c(1,2))
		barplot(snp_remove_by_filter[,1], names.arg = rownames(snp_remove_by_filter), col = terrain.colors(nrow(snp_remove_by_filter), alpha = 1) ,las=2)
		barplot(snp_remove_by_filter[,2], names.arg = rownames(snp_remove_by_filter), col = terrain.colors(nrow(snp_remove_by_filter), alpha = 1) ,las=2, ylim = c(0,1))
		pdf(plot_file_inv)
		par(mar = c(12, 3, 2, 2), cex.axis = 0.6, mfrow = c(1,2))
		barplot(inv_remove_by_filter[,1], names.arg = rownames(inv_remove_by_filter), col = terrain.colors(nrow(inv_remove_by_filter), alpha = 1) ,las=2)
		barplot(inv_remove_by_filter[,2], names.arg = rownames(inv_remove_by_filter), col = terrain.colors(nrow(inv_remove_by_filter), alpha = 1) ,las=2, ylim = c(0,1))
		graphics.off()

		plot_file_snp <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_filters_variable.pdf")
		plot_file_inv <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_filter_invariable.pdf")
		return(list(snp_remove_by_filter = snp_remove_by_filter, inv_remove_by_filter = inv_remove_by_filter, plot_file_snp = plot_file_snp, plot_file_inv = plot_file_inv))
	#plot barplot sites removed
	#add dr.plot for each filter
	}

	create_plot_from_filters <- function (RES, filters, report_dir, chrom.s = 1, chrom.e = NULL, name = "", opt = NULL, ...) {
		if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))
	    create_dr_plot_from_filters <- function(vcf_bed, filter, windows = NULL, pc = TRUE, opt = NULL, ...){
			if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))
	        verbose <- opt$verbose
	        opt$verbose <- FALSE
	        original_vcf <- vcf_bed
	        vcf_bed <- vcf_bed[filter]
	        myDots <- list(...)
	        if (!is.null(myDots$col)) {
	            rbcol <- myDots$col
	            rcol  <- myDots$col
	        }

	        if(sum(filter) == 0){
	            return(list())
	        }else{
	            if(!is.null(windows)){
	                windows_genome <- windows_from_reference(opt$reference, vcf_bed, windows_size = windows)
	                variants_in_windows <- lapply(seq_along(windows_genome), function(x) findOverlaps(vcf_bed, windows_genome[[x]], minoverlap=1, ignore.strand=TRUE))
	                nb_variants_per_windows <- lapply(variants_in_windows, function(y) tapply(y@from, y@to, function(x) length(x) ))[[1]]

	                total_variants_in_windows <- lapply(seq_along(windows_genome), function(x) findOverlaps(original_vcf, windows_genome[[x]], minoverlap=1, ignore.strand=TRUE))
	                total_nb_variants_per_windows <- lapply(total_variants_in_windows, function(y) tapply(y@from, y@to, function(x) length(x) ))[[1]]

	                total_nb_variants_per_windows[names(nb_variants_per_windows)]
	                windows_genome_with_variant <- windows_genome[[1]][as.numeric(names(nb_variants_per_windows))]

	                if(pc) nb_variants_per_windows <- nb_variants_per_windows / total_nb_variants_per_windows[names(nb_variants_per_windows)] * 100
	                values(windows_genome_with_variant) <- data.frame(nb_variant=nb_variants_per_windows)

	                rlst <- as.matrix(data.frame(left=windows_genome_with_variant@ranges@start, bottom=0, right=windows_genome_with_variant@ranges@start + windows_genome_with_variant@ranges@width, top=windows_genome_with_variant@elementMetadata$nb_variant))
	            }else{
	                rlst <- as.matrix(data.frame(left=vcf_bed@ranges@start, bottom=0, right=vcf_bed@ranges@start+vcf_bed@ranges@width, top=1))
	            }
	            drlist <- list(rlst = rlst, chrom.s = 1, chrom.e = chrom.e, title = NULL, rcol = rcol, rbcol = rbcol)
	            return(drlist)
	        }
	    }

	    dr.plot_mod <- function (dmat = NULL, rlst = NULL, chrom.s = 1, chrom.e = NULL,
		    title = NULL, hline = NULL, dcol = NULL, rcol = NULL, rbcol = NULL, opt = NULL,
		    ...){
    			if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))

			    if (class(rlst) != "list" & class(rlst) == "matrix") {
			        rlst <- list(rlst)
			    }
			    if (class(rlst) != "list" & !is.null(rlst)) {
			        stop(paste("parameter rlst is of type", class(rlst),
			            "instead of type list."))
			    }
			    if (is.null(chrom.e)) {
			        stop("chrom.e (end chromosome position) must be specified.")
			    }
			    plot(c(chrom.s, chrom.e), c(0, 0), type = "n", xaxt = "n",
			        xlab = "", ylab = "", las = 1,
			        ...)
			    if (length(rlst) > 1) {
			        rcol <- rep(rcol, times = length(rlst))
			        rbcol <- rep(rbcol, times = length(rlst))
			    }
			    if (!is.null(rlst)) {
			        for (i in 1:length(rlst)) {
			            rmat <- rlst[[i]]
			            graphics::rect(xleft = rmat[, 1], ybottom = rmat[,
			                2], xright = rmat[, 3], ytop = rmat[, 4], col = rcol[i],
			                border = rbcol[i], ...)
			        }
			    }
			    graphics::title(main = title, line = 1)
			    return(invisible(NULL))
		}

	    vcf_bed <- vcf2Grange(vcfRaw(RES), metadata = NULL, opt)
	    if(is.null(chrom.e)) chrom.e = max(vcf_bed@ranges@start)
	    windows_dr_plot  <- 20000
		if(chrom.e < windows_dr_plot * 10) windows_dr_plot <- chrom.e %/% 20

	    mwidth <- 5
	    ncols <- 1
	    nrows <- length(names(filters)) + 1
	    mheight <- 0.5
	    heights <- c()

	    plot_dr_filters_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_", name, "_filters_dr_plot.pdf")

	    pdf(plot_dr_filters_name, width = 35, height = 33)
	    graphics::par(mar=c(3,3,5,3), oma = c(5,5,1,3), cex.axis = 1.5, cex.main=2.3, cex.lab = 2.3 )
	    graphics::layout(matrix(1:c(ncols * nrows), nrow = nrows, ncol = ncols, byrow = TRUE), widths = mwidth, heights = heights)
	    col <- terrain.colors(length(filters), alpha = 1)
	    names(col) <- names(filters)

	    for(i in names(filters)){
	        #print(i)
	        dr_list <- create_dr_plot_from_filters(vcf_bed, filters[[i]], windows = windows_dr_plot, opt = opt, col = col[i])
	        dr.plot_mod(rlst=dr_list$rlst,rcol = dr_list$rcol, rbcol = dr_list$rbcol, chrom.e = chrom.e, title = i, ylim =c(0,110))
	    }

	    dr_list <- create_dr_plot_from_filters(vcf_bed= vcf_bed, filter = rep(TRUE,length(vcf_bed)) , windows = windows_dr_plot, pc = FALSE, opt=opt, col = "darkslategray4" )
	    dr.plot_mod(rlst=dr_list$rlst,rcol = dr_list$rcol, rbcol = dr_list$rbcol, chrom.e = chrom.e, title = "Nucleotides", ylim =c(0,max(unlist(dr_list$rlst[,c(2,4)]))+100))

	    graphics::axis(side = 1, line = 0)
	    graphics::title(xlab = "Base pairs", line = 1.6, outer = TRUE)
	    graphics::title(main = unique(vcf_bed@seqnames), line = -1, outer = TRUE)
	    graphics.off()

	    plot_dr_filters_name <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_", name, "_filters_dr_plot.pdf")
	    opt$verbose <- verbose
	    return(plot_dr_filters_name)
	#plot the filters as bed range
	}

	create_filter_report <- function(opt){
		table_filter_opt <- matrix(NA, ncol = 2, nrow = 12, dimnames = list(c("GATK filter", "Bi-allelic", "Missing", "Depth","Genotype Quality", "Remove indel", "Subtract bedfiles", "Filter by % repeats in windows", "Remove fix heterozygous in all samples","Remove all heterozygous sites", "Remove contiguous fix heterozygous in 50bp windows", "Filter by normalized depth in windows"), c("Active", "Values")))

		table_filter_opt[,1] <- c(
			isTRUE(!is.na(opt$filter_GATK_info)),
			isTRUE(opt$biallelic),
			isTRUE(!is.na(opt$missing)),
			isTRUE(!is.na(opt$filter_depth)),
			isTRUE(!is.na(opt$filter_genotype_quality)),
			isTRUE(opt$filter_indel),
			isTRUE(!is.na(opt$bed_file[1])),
			isTRUE(opt$filter_repeats_by_windows),
			isTRUE(opt$filter_fix_het),
			isTRUE(opt$filter_all_het),
			isTRUE(!is.na(opt$filter_fix_het_contiguous_in_pop)),
			isTRUE(isTRUE(!is.na(opt$filter_high_DP_standardized)))
		)

		High_DP_info <- parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)[5:8]
		DP_info <- unlist(parse_filter_depth(opt$filter_depth)[3:4])

		table_filter_opt[,2] <- c(
		if(isTRUE(is.null(opt$filter_GATK_info))) "" else opt$filter_GATK_info,
		"",
		if(is.na(opt$missing)) "" else sprintf("Proportion of missing genotypes allowed : %s", opt$missing),
		if(is.na(DP_info[1])) "" else sprintf("Setting genotype as missing if min depth < %s and maximum depth >  %s", DP_info["min_DP"], DP_info["max_DP"]),
		if(is.na(opt$filter_genotype_quality[1])) "" else sprintf("Setting minimum genotype quality set as %s for GQ (snps) and %s for RGQ (invariant). ", opt$filter_genotype_quality[["GQ"]], opt$filter_genotype_quality[["RGQ"]]),
		"",
		if(is.na(opt$bed_file[1])) "" else paste(gsub(",", ", ", opt$bed_file), collapse=", "),
		if(!isTRUE(opt$filter_repeats_by_windows)) "" else "Filter 20kb windows with more than 50% bp made up of repeats.",
		"",
		"",
		if(is.na(opt$filter_fix_het_contiguous_in_pop)) "" else sprintf("Remove 50bp windows with at least 2 contiguous fixed heterozygous in the pop matching the pattern %s", opt$filter_fix_het_contiguous_in_pop),
		if (is.na(High_DP_info)[1]) "" else sprintf("Remove %s bp windows with an step of %s bp, which have a median standardized depth > %s and a 95 quartile > %s across all samples",High_DP_info$windows_size, High_DP_info$slidding, High_DP_info$threshold[1], High_DP_info$threshold_CI )
		)

		return(table_filter_opt)
	}

	#plot indel size
	plot_indel_size <- function(RES, report_dir, vcf_file){
		vcf_indel <- vcfRaw(RES)[Filters(RES)$indel]
		get_indel_size <- function(insertion){
			insertion[sapply(insertion, nchar) == 1 & insertion != "*"] <- NA
			a <- strsplit(insertion, ",")
			b <- sapply(a, function(x) max(nchar(x)))
			insertion[b == 1 & !grepl("\\*", insertion)] <- NA
			indel_size <- b[!is.na(insertion)]
			return(indel_size)
		}

		insertion <- get_indel_size(getALT(vcf_indel))
		deletion  <- get_indel_size(getREF(vcf_indel))
		indel_size_plot_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(basename(vcf_file))), "_indel_size.pdf")
		if(nrow(vcf_indel) != 0){
			pdf(indel_size_plot_name)
			hist(c(insertion, deletion), 100, xlab = "", col ="aquamarine2", main ="")
			graphics.off()
			indel_size_plot_name <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(basename(vcf_file))), "_indel_size.pdf")
		}
		return(indel_size_plot_name)
	}
	#recreate the plots of filter vs non filter stats if the file were splitted

	#Plot chromose
	create_chromqc_plots <- function(RES, chr = unique(getCHROM(vcfRaw(RES))), report_dir = ".", opt = NULL){
		if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))
		void <- graphics.off()
		dna <- ape::as.DNAbin(reference(RES)[names(reference(RES)) == chr])
		chrom <- create.chromR(name=chr, vcf=vcfRaw(RES)[getCHROM(vcfRaw(RES)) == chr], seq=dna, verbose=FALSE)
		chrom_filtered <- create.chromR(name=chr, vcf=vcf(RES)[getCHROM(vcf(RES)) == chr], seq=dna, verbose=FALSE)

		chrom <- proc.chromR(chrom, win.size = 20000, verbose = FALSE)
		chrom_filtered <- proc.chromR(chrom_filtered, win.size = 20000, verbose = FALSE)

		# plot_chrom_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "chrom_plot.pdf")
		# plot_chrom_name_filt <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "chrom_filtered_plot.pdf")
		# pdf(plot_chrom_name)
		# plot(chrom)
		# pdf(plot_chrom_name_filt)
		# plot(chrom_filtered)
		# graphics.off()

		plot_chromqc_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)),"_", chr, "_chrom_plot_qc.pdf")
		plot_chromqc_name_filt <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz","",(vcf_name)),"_", chr, "_chrom_filtered_plot_qc.pdf")
		pdf(plot_chromqc_name, compress = TRUE)
		try(chromoqc(chrom, xlim =c(1, length(dna[[1]]))))
		graphics.off()
		pdf(plot_chromqc_name_filt, compress = TRUE)
		chromoqc(chrom_filtered, xlim =c(1, length(dna[[1]])))
		graphics.off()

		# plot_chrom_name <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "chrom_plot.pdf")
		#plot_chrom_name_filt <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)),"_", chr, "_chrom_filtered_plot.pdf")
		plot_chromqc_name <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_", chr, "_chrom_plot_qc.pdf")
		plot_chromqc_name_filt <- paste0("./plots/", gsub(".vcf|.vcf.gz","",(vcf_name)), "_", chr, "_chrom_filtered_plot_qc.pdf")
		return(list(plot_chromqc_name_filt=plot_chromqc_name_filt, plot_chromqc_name=plot_chromqc_name ))
	}

	create_plot_SFS <- function(RES, report_dir, vcf_name, pop = NA){
		gt <- extract.gt(vcf(RES), convertNA = FALSE)
		nas 	<- rowSums(gt == "./." | gt == ".|.")
		count_het 	<- rowSums(gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0" , na.rm = TRUE)
		count_hom_alt 	<- rowSums(gt == "1/1" | gt == "1|1" , na.rm = TRUE)
		sum_ALT_allele <- count_het + count_hom_alt * 2
		#we need to calculate the frequency
		#considering the missing data.
		# We will count the number of alternate alleles
		# and divide by the number of samples -numbers of NA
		# Then we round it and multiply again by the number of sample

		x <- sum_ALT_allele / ((ncol(gt) - nas)*2)
		# find the closest value of possible frequency : 0 to max nb of individual
		possible_freq <- seq(0, ncol(gt)*2, by = 1) / (ncol(gt)*2)
		all_freq <- unique(x)

		#Create a lookup table to assign frequency with NA to the right frequency
		lookup <- data.frame(find= all_freq, replace = NA)
		for(i in 1:length(all_freq)){
			lookup[i,"replace"] <- possible_freq[which.min(abs(possible_freq -  all_freq[i]))]
		}
		for(i in 1:nrow(lookup)){
			x[x == lookup[i, "find"]] <- lookup[i, "replace"]
		}
		sum_ALT_allele <- x * (ncol(gt)*2)
		#table of contengency counting each site in each frequency(), we are using tabulate instead of table to set the correct number of bins even if they are empty. Here the SFS must have  numbers of chromosome minus missing (they all have now the same number of chromosomes after the random adding of NA) and 1 extra bin for non-variant.
		SFS <- tabulate(as.integer(round(sum_ALT_allele)), ncol(gt)*2)
		#add no variation sites (a nb of sites)
		no_var <- sum(width(invariant(RES)))
		SFS <- c(no_var, SFS)

		#add the fixed alternates sites
		fixed_ALT <- sum(Filters(RES)$inv_fix_sites)
		SFS[length(SFS)] <- fixed_ALT

		fSFS <- (SFS +rev(SFS))[1:ceiling(length(SFS)/2)]
		fSFS[ceiling(length(SFS)/2)] <- fSFS[ceiling(length(SFS)/2)] /2

		down_SFS <- downsample_SFS(SFS, 10)
		down_fSFS <- downsample_SFS(fSFS, 10)

		make_xlabels <- function(SFS, by  = 3, proportion = FALSE){
			at <- seq(1, length(SFS)-1, by =by) - 0.5
			if(proportion){
				xlabels <- round(seq(1, length(SFS)-1, by  = by) / (length(SFS)-1), 2)
			}else{
				xlabels <- seq(1, length(SFS)-1, by  = by)
			}
				return(list(at=at, xlabels=xlabels))
		}
		make_second_yaxis <- function(SFS){
			at <- seq(0 , max(SFS[-1]), by =  max(SFS[-1]) %/%10 )
			prop_SFS <- SFS[-1]/sum(SFS[-1])
			YLIM <- round(seq(0, max(prop_SFS), by = max(prop_SFS) / 10 ), 2)
			if(length(at) != length(YLIM)){
				YLIM[length(at)] <- max(YLIM) + round(max(prop_SFS) / 10,2)
			}
			return(list(at=at, YLIM=YLIM))
		}

		plot_SFS <- function(SFS, by = 2, ...){
			barplot(SFS[-1] ,space = 0, las = 1, ylim = c(0, max(SFS[-1])), xaxt="n", ...)
			axis(1, line = 0 , at = make_xlabels(SFS, by)$at, tick = TRUE, labels = make_xlabels(SFS, by)$xlabels)
			axis(1, line = 2 , at = make_xlabels(SFS, by, proportion = TRUE)$at, tick = TRUE, labels = make_xlabels(SFS, by, proportion = TRUE)$xlabels)
			try(axis(4, line = -1 , at = make_second_yaxis(SFS)$at, tick = TRUE, labels = make_second_yaxis(SFS)$YLIM))
		}

		plot_SFS_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz", "",(vcf_name)), "_SFS.pdf")
		plot_downSFS_name <- paste0(report_dir, "/plots/", gsub(".vcf|.vcf.gz", "",(vcf_name)), "_SFS_donsampled10.pdf")

		pdf(plot_SFS_name, width = 50)
		par(mfrow = c(2,1))
		plot_SFS(SFS, main = "Unfolded SFS", col = "aquamarine2")
		plot_SFS(fSFS, main = "Folded SFS", col = "aquamarine2")
		graphics.off()

		pdf(plot_downSFS_name, 25)
		par(mfrow = c(2,1))
		tryCatch(plot_SFS(down_SFS, by = 1, main = "Unfolded SFS downsampled to 10" , col = "aquamarine2"), error= function(e) message(sprintf("\n\nPlot for downsampled unfolded SFS failed, downsample chromosome nb ,%i, might be lower than actual number of unfolded chromosomes", 10)))
		tryCatch(plot_SFS(down_fSFS, by = 1, main = "Folded SFS downsampled to 10" , col = "aquamarine2"), error= function(e) message(sprintf("\n\nPlot for downsampled folded SFS failed, downsample chromosome nb ,%i, might be lower than actual number of folded chromosomes (half)", 10)))
		graphics.off()

		plot_SFS_name <- paste0("./plots/", gsub(".vcf|.vcf.gz", "",(vcf_name)), "_SFS.pdf")
		plot_downSFS_name <- paste0("./plots/", gsub(".vcf|.vcf.gz", "",(vcf_name)), "_SFS_donsampled10.pdf")

		mat_SFS <- matrix("", ncol=4, nrow=length(SFS), dimnames = list(0:(length(SFS)-1), c("Unfolded SFS", "Folded SFS", "Downsampled unfolded SFS", "Downsampled folded SFS")))
		mat_SFS[,1] <- SFS
		mat_SFS[1:length(fSFS),2] <- fSFS
		mat_SFS[1:length(down_SFS),3] <- down_SFS
		mat_SFS[1:length(down_fSFS),4] <- down_fSFS

		write.table(t(mat_SFS), paste0(report_dir, "/",  gsub(".vcf|.vcf.gz", "",(vcf_name)), ".sfs"))

		return(list(SFS=SFS, fSFS=fSFS,down_SFS=down_SFS, down_fSFS=down_fSFS, plot_SFS_name=plot_SFS_name, plot_downSFS_name=plot_downSFS_name, mat_SFS = as.data.frame(mat_SFS)))
	}

	RES_by_chr <- function(RES, chr=NULL){
		if(is.null(chr)){
			chr <- unique(getCHROM(vcfRaw(RES)))
		}

		res_vcf <- as.list(rep(NA, length(chr)))
		names(res_vcf) <- chr

		for(i in chr){
			index_pos_chrom <- getCHROM(vcfRaw(RES)) == i
			res_vcf[[i]]$vcf <- suppressWarnings(vcfRaw(RES)[index_pos_chrom,])
			res_vcf[[i]]$master_filter <- masterFilter(RES)[index_pos_chrom]
			res_vcf[[i]]$filters <- lapply(Filters(RES), function(x) return(x[index_pos_chrom]))
			res_vcf[[i]]$filter_DP_repeats[[1]] <- lapply(Normalized_DP_Filter(RES)[[1]], function(x) return(x[as.logical(x@seqnames == i)]))
			res_vcf[[i]]$filter_DP_repeats[[2]] <- lapply(Normalized_DP_Filter(RES)[[2]], function(x) return(x[as.logical(x@seqnames == i)]))
		}
		return(res_vcf)
	}

	matrix_filter_intersection <- function(filters) {
		Filters <- filters[sapply(filters, class) == "logical"]
		mat_and <- matrix(NA, ncol = length(Filters), nrow = length(Filters))
		mat_or <- matrix(NA, ncol = length(Filters), nrow = length(Filters))
		for(i in seq_along(Filters)){
			for(j in seq_along(Filters)){
				mat_and[i,j] <- sum(Filters[[i]] & Filters[[j]])
				mat_or[i,j] <- sum(Filters[[i]] | Filters[[j]])
			}
		}

		dimnames(mat_and) <- dimnames(mat_or) <- list(names(Filters), names(Filters))
		prop_mat_and <- round(mat_and / length(Filters[[1]]), 2)
		prop_mat_or <- round(mat_or / length(Filters[[1]]), 2)
		return(list(filter_mat_and = mat_and, filter_prop_mat_and = prop_mat_and, filter_mat_or = mat_or, filter_prop_mat_or = prop_mat_or))
	}

	Create_HTML_report = function(file = "report.html", append = TRUE, directory = getwd(), opt = NULL) {
		if(is.null(opt)) opt <- parse_args(OptionParser(option_list=initialise_option()))
		wd <- getwd()
		setwd(directory)
		file = file.path(directory, file)
		if(file.exists(file)) file.remove(file)

		HTMLInitFile(outdir = directory, filename=basename(gsub(".html", "",file)), extension="html",
		HTMLframe=FALSE, BackGroundColor = "FFFFFF", BackGroundImg = "#A0A0A0",
		Title = "", CSSFile=system.file("samples", "R2HTML.css", package="R2HTML"), useLaTeX=TRUE, useGrid=TRUE)

		file <- paste0(file, ".html")
		#cat("\n", file = file, append = append)
		HTML.title("Reports VCF filtering",file = file, HR=1, align = "center")
		HTML(sprintf("VCFfile : %s", vcf_name),file = file)
		HTML(sprintf("Reference : %s", basename(opt$reference)),file = file)
		HTML(opt$call,file = file)

		#add call of function with arguments

		#info on samples
		HTMLhr(file = file, Width = "100%", Size = "2")
		HTML(suppressWarnings(matrix(colnames(vcf(RES)@gt)[-1], ncol = 4 )),file = file, Border = 1, caption = sprintf("List of samples (%s)", ncol(vcf(RES)@gt)-1), captionalign = "top")
		HTMLhr(file = file)

		#general info on filtering
		HTML.title("Global filtering",file = file, HR=1)
		HTML("Table of sites filtered\n NOTE: variable sites can still contain indels depending if the indel filter has been applied or not.", file = file)
		HTML(opt$nb_sites_reports$table_report, file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTMLInsertGraph(opt$nb_sites_reports$plot_file , Caption = "Summary of filtered sites", file = file, GraphBorder=1, Align="left", WidthHTML=500, HeightHTML=NULL)

		HTMLhr(file = file)
		HTML.title("Info on quality of the filtered vcf",file = file, HR=1)
		HTMLhr(file = file)

		# add multiple plot if multiple chrom !!!
		if(isTRUE(opt$single_chrom)){
			HTMLInsertGraph(opt$chromqc$plot_chromqc_name , Caption = "Original", file = file, GraphBorder=1, Align="left", WidthHTML=1500, HeightHTML=NULL)
			HTMLInsertGraph(opt$chromqc$plot_chromqc_name_filt , Caption = "Filtered", file = file, GraphBorder=1, Align="left", WidthHTML=1500, HeightHTML=NULL)
		}else{
			for(i in unique(getCHROM(vcfRaw(RES)))){
				HTMLInsertGraph(opt$chromqc[[i]]$plot_chromqc_name , Caption = "Original", file = file, GraphBorder=1, Align="left", WidthHTML=1500, HeightHTML=NULL)
				HTMLInsertGraph(opt$chromqc[[i]]$plot_chromqc_name_filt , Caption = "Filtered", file = file, GraphBorder=1, Align="left", WidthHTML=1500, HeightHTML=NULL)
			}
		}

		#info on single filter impact
		HTMLhr(file = file)
		HTML.title("Filters applied",file = file, HR=1)
		HTML("",file = file)

		HTML(create_filter_report(opt) ,file = file, border =2,innerBorder = 2, caption = "List of the different filters applied", captionalign = "top")

		plot_filter <- create_plot_filter(RES, report_dir, opt)
		HTML(plot_filter$snp_remove_by_filter,file = file, Border = 1, caption = "Summary of filtered variable sites by each filter", captionalign = "top")
		HTMLInsertGraph(plot_filter$plot_file_snp , Caption = "Summary of filtered variable sites by each filter", file = file, GraphBorder=1, Align="center", WidthHTML=500, HeightHTML=NULL)
		HTMLhr(file = file)
		HTML(plot_filter$inv_remove_by_filter,file = file, Border = 1, caption = "Summary of filtered invariable sites by each filter", captionalign = "top")
		HTMLInsertGraph(plot_filter$plot_file_inv , Caption = "Summary of filtered invariable sites by each filter", file = file, GraphBorder=1, Align="center", WidthHTML=500, HeightHTML=NULL)

		HTMLhr(file = file)
		filter_matrix_snp <- matrix_filter_intersection(Filters(RES))
		filter_matrix_inv <- matrix_filter_intersection(FiltersInv(RES))
		HTML("Tables for SNP",file = file)
		HTML("Tables showing the number and percentage of sites that are common in two filters",file = file)
		HTML((filter_matrix_snp$filter_mat_and), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML((filter_matrix_snp$filter_prop_mat_and), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML("Tables showing the number and percentage of sites that are additive in two filters",file = file)
		HTML((filter_matrix_snp$filter_mat_or), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML((filter_matrix_snp$filter_prop_mat_or), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")

		HTML("Tables for invariant",file = file)
		HTML("Tables showing the number and percentage of sites that are common in two filters",file = file)
		HTML((filter_matrix_inv$filter_mat_and), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML((filter_matrix_inv$filter_prop_mat_and), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML("Tables showing the number and percentage of sites that are additive in two filters",file = file)
		HTML((filter_matrix_inv$filter_mat_or), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")
		HTML((filter_matrix_inv$filter_prop_mat_or), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")

		#plots of the sites removed by each filter
		if(isTRUE(opt$single_chrom)){
			HTMLInsertGraph(opt$dr_plot_filter, Caption = "Percentage of sites removed by each filter in 20kb windows", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
		}else{
			for( i in unique(getCHROM(vcfRaw(RES)))){
				HTMLInsertGraph(opt$dr_plot_filter[[i]], Caption = "Percentage of sites removed by each filter in 20kb windows", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
			}
		}
		HTMLhr(file = file)
		HTMLhr(file = file)

		if(isTRUE(!is.na(opt$filter_high_DP_standardized))){
			if(isTRUE(opt$single_chrom)){
				HTML.title("Impact of different thresholds for high normalized depth by windows", file = file, HR=3)
				HTML("These filters are designed to try to deal with repetitive genome. Region with repeats will often have higher depth due to read mapping to two or more regions. We use short windows (recommended between 1000 and 5000) to calculate the depth for each samples for each snps in the VCFfile. The depth is nomralized for each sample by the median depth over all snps. This enable to compare samples with different coverage and apply relative threshold instead of raw value of depth. Then two kind of rules are applied.\nFirst a windows is removed if the mean normalized depth is higher than the threshold and the mean upper 95 quartile is above the CI threshold. For example if the a sample has median depth of 30X, a windows with 60X and a 95 quartile of 90X will have a normalized median depth of 2 and 95 quartile of 3. Taking the mean over samples for this windows might give an value of normalized median depth of 2.2 and 95 quartile of 3.2. We than apply the threshold (typically 2 and 4) 2.2 > 2 or 3.2 > 4, hence the snps in this windows are filtered.\nThe second rule is more stringent and remove a windows if any sample is over the thresholds for a windows.\nThe pipeline use the first rule with the first value specified in the option --filter_high_DP_standardized",file = file)
				HTMLInsertGraph(opt$dr_plot_filter_high_DP_across , Caption = "Percentage of sites removed by each normalized depth filter using the first rule", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
				HTMLInsertGraph(opt$dr_plot_filter_high_DP_total , Caption = "Percentage of sites removed by each normalized depth filter using the second rule", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
			}else{
				void <- lapply(opt$dr_plot_filter_high_DP_across, HTMLInsertGraph, Caption = "Percentage of sites removed by each normalized depth filter using the first rule", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
				void <- lapply(opt$dr_plot_filter_high_DP_total, HTMLInsertGraph, Caption = "Percentage of sites removed by each normalized depth filter using the second rule", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
			}
		}

		HTMLbr(x=1,file =file,append=TRUE)
		HTMLhr(file = file)
		HTML.title("Info on indels",file = file, HR=1)
		HTMLInsertGraph(opt$indel_plots , Caption = "Size of indels", file = file, GraphBorder=1, Align="left", WidthHTML=500, HeightHTML=NULL)

		HTMLhr(file = file)
		#Add SFS
		HTML.title("Site frequency spectrum after the filtering",file = file, HR=1)
		HTML("Note that missing data are taken into account. The approximation is correct if there is not too much missing allowed (~20%) controled by the option --missing", file = file)
		HTMLInsertGraph(opt$SFS_report$plot_SFS_name , Caption = "", file = file, GraphBorder=1, Align="center", WidthHTML=1500, HeightHTML=NULL)
		HTMLInsertGraph(opt$SFS_report$plot_downSFS_name , Caption = "", file = file, GraphBorder=1, Align="center", WidthHTML=500, HeightHTML=NULL)


		HTML(t(opt$SFS_report$mat_SFS), file = file, Border = 1, innerBorder = 1, classfirstline = "header",
		   classfirstcolumn = "header", append = TRUE, align = "left", caption = "", captionalign = "bottom",
		   classcaption = "captiondataframe", classtable = "dataframe")

		HTMLhr(file = file)
		HTMLhr(file = file)
		HTML("Program written by Benjamin Laenen", file = file)

		r_logo <- "https://upload.wikimedia.org/wikipedia/commons/1/1b/R_logo.svg"
		if (!file.exists("R_logo.svg")) suppressWarnings(try(download.file(r_logo, "R_logo.svg", quiet = TRUE)))
		void <- HTMLInsertGraph("R_logo.svg", file = file, Align = "left", WidthHTML=50)

		message(paste("Report written: ", file, sep = ""), appendLF=TRUE)
		setwd(wd)
	}

	report_dir <- paste(outputdir, "report", sep = "/")
	report_plots <- paste(report_dir, "plots", sep = "/")
	if(!dir.exists(outputdir)) dir.create(outputdir)
	if(!dir.exists(report_dir)) dir.create(report_dir)
	if(!dir.exists(report_plots)) dir.create(report_plots)

	opt$indel_plots <- plot_indel_size(RES, report_dir, vcf_file=opt$vcf_file)

	#loop over chr if mutiple
	if(isTRUE(opt$single_chrom)){
		if(isTRUE(skip_chromqc)){
			opt$chromqc <- list(plot_chromqc_name_filt = NA, plot_chromqc_name = NA)
		}else{
			opt$chromqc <- create_chromqc_plots(RES, unique(getCHROM(vcfRaw(RES))), report_dir, opt = opt)
			#the chromqc are a bit large and slow down the process
			opt$chromqc <- lapply(opt$chromqc, function(x) gsub(".pdf", ".No.pdf", x))
		}

		opt$dr_plot_filter <- create_plot_from_filters(RES, Filters(RES), report_dir=report_dir, opt = opt)
		
		if(all(!is.na(opt$filter_high_DP_standardized))){
			# Plot the result from all normalized DP filters
			# add if
			vcf_bed_GRange <- vcf2Grange(vcfRaw(RES))
			names(Normalized_DP_Filter(RES)[[1]]) <- names(Normalized_DP_Filter(RES)[[2]]) <- paste0("median_normDP_", parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)$threshold)
			filters_high_DP_standardized_across<- lapply(Normalized_DP_Filter(RES)[[2]], function(x) Grange2filter(x, vcf_bed_GRange,opt))
			filters_high_DP_standardized_total <- lapply(Normalized_DP_Filter(RES)[[1]], function(x) Grange2filter(x, vcf_bed_GRange,opt))

			opt$dr_plot_filter_high_DP_across <- create_plot_from_filters(RES, filters_high_DP_standardized_across, name = "_high_DP_across",  report_dir=report_dir, opt = opt)
			opt$dr_plot_filter_high_DP_total <- create_plot_from_filters(RES, filters_high_DP_standardized_total, name = "_high_DP_per_sample", report_dir=report_dir, opt = opt)
		}
	}else{
		for(i in unique(getCHROM(vcfRaw(RES)))){
			#Loop over chr
			opt$chromqc[[i]] <- create_chromqc_plots(RES, i, report_dir, opt = opt)
			#file.rename(unlist(opt$chromqc[[i]]), gsub(".pdf", paste0("_", i, ".pdf"), opt$chromqc[[i]]))
			#opt$chromqc[[i]] <- lapply(opt$chromqc[[i]], function(x) gsub(".pdf", paste0("_", i, ".pdf"), x))

			VCF_chr <- suppressWarnings(RES_by_chr(RES, chr = i))
			opt$dr_plot_filter[[i]] <- create_plot_from_filters(VCF_chr[[i]], VCF_chr[[i]]$filters, report_dir=report_dir, name = i)

			if(isTRUE(!is.na(opt$filter_high_DP_standardized))){
				VCF_chr_bed_GRange <- vcf2Grange(VCF_chr[[i]]$vcf)
				names(VCF_chr[[i]]$filter_DP_repeats[[1]]) <- names(VCF_chr[[i]]$filter_DP_repeats[[2]]) <- paste0("median_normDP_", parse_filter_high_DP_standardized(opt$filter_high_DP_standardized)$threshold)
				filters_high_DP_standardized_across<- lapply(VCF_chr[[i]]$filter_DP_repeats[[1]], function(x) Grange2filter(x, VCF_chr_bed_GRange,opt))
				filters_high_DP_standardized_total <- lapply(VCF_chr[[i]]$filter_DP_repeats[[2]], function(x) Grange2filter(x, VCF_chr_bed_GRange,opt))

				opt$dr_plot_filter_high_DP_across[[i]] <- create_plot_from_filters(VCF_chr[[i]], filters_high_DP_standardized_across, name = paste0("_", i, "_high_DP_across"),  report_dir=report_dir)
				opt$dr_plot_filter_high_DP_total[[i]] <- create_plot_from_filters(VCF_chr[[i]], filters_high_DP_standardized_total, name = paste0("_", i, "_high_DP_per_sample"), report_dir=report_dir)
			}
		}
	}

	opt$nb_sites_reports <- make_plot_nb_sites(RES, report_dir, opt)
	opt$SFS_report <- create_plot_SFS(RES, report_dir, vcf_name)
	opt$call <- gsub(" -", " \\\\\n-" ,get_call())
	out <- Create_HTML_report(file = opt$report, append = TRUE, directory = normalizePath(report_dir), opt = opt)
}
