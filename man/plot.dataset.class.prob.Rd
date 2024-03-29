% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotClass.R
\name{plot.dataset.class.prob}
\alias{plot.dataset.class.prob}
\title{Plots a stacked bar chart of most likely latent class proportion, given a matrix of covariates.}
\usage{
\method{plot}{dataset.class.prob}(X, alpha, title = "Proportion of Latent Classes")
}
\arguments{
\item{X}{A matrix of covariates.}

\item{alpha}{A matrix of logit regression coeffiecients.}

\item{title}{A text string for plot title.}
}
\value{
A \code{ggplot2} object.
}
\description{
Plots a stacked bar chart of most likely latent class proportion, given a matrix of covariates.
}
\seealso{
\code{\link{LRMoE.fit}}, \code{\link{predict.class}}
}
