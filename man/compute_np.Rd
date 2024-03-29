% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/np_pipeline.R
\name{compute_np}
\alias{compute_np}
\title{main function to compute np from a user-provided expression matrix.}
\usage{
compute_np(
  cache = NULL,
  experiment_name,
  ppi = "biogrid",
  min_score = NULL,
  exp_mat,
  mir_paper = TRUE,
  ncores = 1
)
}
\arguments{
\item{cache}{user-provided filepath for where to store data etc}

\item{experiment_name}{name of the experiment for saving output.}

\item{ppi}{should we use biogrid or stringdb for the PPI}

\item{min_score}{if ppi is stringdb, which mininum score should we use to filter edges?}

\item{exp_mat}{expression matrix where columns are samples and rows are features}

\item{mir_paper}{are we running this in the context of the mir paper? a few quirks of that data}

\item{ncores}{number of cores to use for calculations}
}
\value{
tidy data frame with one column for expression and another for np
}
\description{
main function to compute np from a user-provided expression matrix.
}
