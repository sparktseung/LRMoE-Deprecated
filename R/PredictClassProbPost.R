#' Predict the posterior latent class probabilities, given a fixed covariate matrix X, Y and a model.
#'
#' @param X A matrix of covariates.
#' @param Y A matrix of observed responses for \code{X}.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of latent class probabilities by observation and by component.
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @export predict.class.prob.posterior
predict.class.prob.posterior = function(X, Y, alpha, comp.dist, zero.prob, params.list)
{
  gate.ll = gate.logit(X, alpha)
  expert.list = expert.loglik.dim.comp(Y, comp.dist, zero.prob, params.list)
  ll.list = gate.expert.loglik(alpha, gate.ll, expert.list, penalty=FALSE) # , hyper.alpha, hyper.params)
  weighting = comp.zkz.e.recur(gate.ll, expert.list, ll.list)$z.e.obs

  return( weighting )
}
