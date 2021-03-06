% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disease_enrichment.R
\name{gene_to_NCBI}
\alias{gene_to_NCBI}
\title{gene to NCBI}
\usage{
gene_to_NCBI(
  genes,
  organism = "mouse",
  update_info = FALSE,
  parallel = NULL,
  conf_file = getOption("RMINC_BATCH_CONF")
)
}
\arguments{
\item{genes}{gene acronyms}

\item{parallel}{how many processors to run on (default=single processor). Specified as a two element vector, with the first element corresponding to the type of parallelization, and the second to the number of processors to use. For local running set the first element to "local" or "snowfall" for back-compatibility, anything else will be run with batchtools see pMincApply Leaving this argument NULL runs sequentially.}

\item{conf_file}{A batchtools configuration file defaulting to getOption("RMINC_BATCH_CONF")}
}
\value{
NCBI gene IDs
}
\description{
Converts gene acronym to NCBI gene IDs
}
