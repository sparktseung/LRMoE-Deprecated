% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatorSet.R
\name{dataset.simulator}
\alias{dataset.simulator}
\title{Simulate y, given a fixed covariate matrix X and a model}
\usage{
dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)
}
\arguments{
\item{X}{A N*P matrix of of covariates. The first column must be 1. Each row may be different.}

\item{alpha}{A g*P matrix of logit regression coefficients.}

\item{comp.dist}{A d*g matrix of strings, describing component distributions by dimension and by component.}

\item{zero.prob}{A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.}

\item{params.list}{A list of length d, where each element is a sublist of length g.
Each sublist contains one numeric vector, which is the initial parameter guess for the corresponding \code{comp.dist}.}
}
\value{
A matrix of simulated values, where each row represents a policyholder and each column a dimension of the response variable.
}
\description{
Simulate y, given a fixed covariate matrix X and a model
}
