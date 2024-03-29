% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EMalpha.R
\name{alpha.m.recur}
\alias{alpha.m.recur}
\title{ECM: M-Step for logit regression coefficients \code{alpha}.}
\usage{
alpha.m.recur(X, alpha, comp.zkz.e.list, alpha.iter.max, penalty, hyper.alpha)
}
\arguments{
\item{X}{A N*P matrix of numerical covariates.}

\item{alpha}{A g*P matrix of old logit regression coefficients.}

\item{comp.zkz.e.list}{An object returned by \code{\link{comp.zkz.e.recur}}.}

\item{alpha.iter.max}{Numeric: maximum number of iterations.}

\item{penalty}{TRUE/FALSE, which indicates whether penalty is applied.}

\item{hyper.alpha}{A numeric of penalty applied to \code{alpha}.}
}
\value{
\code{alpha.new} Updated logit regression coefficients.
}
\description{
ECM: M-Step for logit regression coefficients \code{alpha}.
}
\keyword{internal}
