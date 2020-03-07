#' Predict the most likely latent class, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A vector of the most likely latent class by observation.
#'
#' @rawNamespace S3method(predict, class)
#'
#' @export predict.class
predict.class = function(X, alpha)
{
  X.alpha = X %*% t(alpha)
  result = matrix(apply(X.alpha, 1, FUN = "which.max"))
  return( result )
}
