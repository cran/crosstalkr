% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_np.R
\name{calc_np_all}
\alias{calc_np_all}
\title{function to calculate the network potential for each protein in a user-provided vector - cpp internal version}
\usage{
calc_np_all(exp, g, v = "default", neighbors = NULL)
}
\arguments{
\item{exp}{expression vector - assumed to be a named vector where the values are expression and the names are the gene name}

\item{g}{igraph object - will be filtered so that only nodes found in both exp and g are kept}

\item{v}{character vector of nodes over which to calculate network potential.}

\item{neighbors}{named list containing the neighbors for each node of graph g. If not provided,
it will be computed}
}
\value{
dataframe containing network potential for each of the inputed gene names.
}
\description{
function to calculate the network potential for each protein in a user-provided vector - cpp internal version
}
