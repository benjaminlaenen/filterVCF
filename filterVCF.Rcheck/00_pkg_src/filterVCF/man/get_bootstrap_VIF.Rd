% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVCF_vif_functions.R
\name{get_bootstrap_VIF}
\alias{get_bootstrap_VIF}
\title{Get one bootstrap replicates of VIF}
\usage{
get_bootstrap_VIF(vif, index = 1)
}
\arguments{
\item{vif}{a VIF object}

\item{index}{index of a single bootstrap to extract, should be lower than NBoot(vif)}
}
\value{
a new VIF object that has been bootstrapped. Note that is not a
proper VIF and will not be compatible with all functions (e.g. intersect)
because it contains duplicated lines as it is a bootstrap.
}
\description{
Bootstrap is not performed directly when using bootstrap_VIF() to avoid
overflowing the memmory with repeated objects. Instead index to create
bootstrap repliactes are created in order to get boostraped VIF on the fly.
This function retrieved one boostrap using the index created in the first
step.
}
\examples{

#create the index for 20kb windows bootstrap
vif <- bootstrap_VIF(vif, nb_bootstrap = 200, bootstrap_by = 20000)
#get the first replicate
boot_vif <- list()
boot_vif[[1]]  get_bootstrap_VIF(vif, index = 1)
#extract SFS from boostrap replicate
get_SFS_vif(boot_vif[[1]])

}
\seealso{
objects to See Also as \code{\link{help}},
}
\author{
~~Benjamin Laenen~~
}
\keyword{~bootstrap}
