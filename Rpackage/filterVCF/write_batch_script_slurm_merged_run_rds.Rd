\name{write_batch_script_slurm_merged_run_rds}
\alias{write_batch_script_slurm_merged_run_rds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{write_batch_script_slurm_merged_run_rds
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write_batch_script_slurm_merged_run_rds(vcf_file, outputdir, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf_file}{
%%     ~~Describe \code{vcf_file} here~~
}
  \item{outputdir}{
%%     ~~Describe \code{outputdir} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vcf_file, outputdir, ...)
{
    bin_bash <- sprintf("#!/bin/bash -l\n#SBATCH -A snic2017-7-175\n#SBATCH -J filter_vcf\n#SBATCH -p core -n \%s\n#SBATCH -t \%s\n\ncd \%s\n",
        opt$no_cores, "60:00:00", outputdir)
    cmd <- get_call()
    cmd <- gsub("--vcf_file .+vcf??? -|-I .+ -", paste0("--vcf_file ",
        basename(vcf_file), " -"), cmd)
    cmd_keep_rds_only <- ""
    cmd <- gsub("--dry_run", "", cmd)
    cmd <- paste(cmd, " --merge_only")
    bash_job_file <- sprintf("merge_run_rds_\%s.sh", gsub(".vcf.gz|.vcf",
        "", basename(vcf_file)))
    void <- write(c(bin_bash, cmd, cmd_keep_rds_only), bash_job_file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
