% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GammaCountNew.R
\name{mgammacount}
\alias{mgammacount}
\title{Calculates moments of \code{\link[rmutil]{GammaCount}} using finite approximation.}
\usage{
mgammacount(order, m, s, tol = 1e-10)
}
\arguments{
\item{order}{A vector of positive power indices.}

\item{m, s}{Paramaters of Gamma count distribution.}

\item{tol}{Cut-off probability threshold. Values above (1-tol) are discarded.}
}
\value{
A vector of Gamma count distribution moments.
}
\description{
Calculates moments of \code{\link[rmutil]{GammaCount}} using finite approximation.
}
\seealso{
\code{\link[rmutil]{GammaCount}}.
}
