#' Predict the most likely latent class, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param Y A matrix of observed responses for \code{X}.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.init A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A vector of the most likely latent class by observation.
#'
#' @export predict.class.posterior
predict.class.posterior = function(X, Y, alpha, comp.dist, zero.prob, params.list)
{
  X.alpha = predict.class.prob.posterior(X, Y, alpha, comp.dist, zero.prob, params.list)
  result = matrix(apply(X.alpha, 1, FUN = "which.max"))
  return( result )
}
