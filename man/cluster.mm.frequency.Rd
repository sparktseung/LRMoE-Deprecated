% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitParams.R
\name{cluster.mm.frequency}
\alias{cluster.mm.frequency}
\title{Initializes parameter for frequency distributions using CMM.}
\usage{
cluster.mm.frequency(Y, cluster)
}
\arguments{
\item{Y}{A vector of response variables.}

\item{cluster}{The \code{cluster} list vector returned by \code{\link[stats]{kmeans}}}
}
\value{
A list of parameter initialization.
}
\description{
Initializes parameter for frequency distributions using CMM.
}
