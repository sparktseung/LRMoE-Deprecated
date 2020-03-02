#' Predict the latent class probabilities, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of latent class probabilities by observation and by component.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @export predict.class.prob
predict.class.prob = function(X, alpha)
{
  weighting = exp(gate.logit(X, alpha))
  return( weighting )
}
