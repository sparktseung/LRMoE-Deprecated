#' A function to calculate the limited mean of a matrix of given distributions (positive part only), by dimension and by component.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param limit A vector of limit values by dimension.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of limited mean values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export dim.comp.limit.y.pos
dim.comp.limit.y.pos = function(comp.dist, params.list, limit)
{
  n.comp = ncol(comp.dist)
  dim.m = nrow(comp.dist)

  result = matrix(0, nrow = dim.m, ncol = n.comp)

  for(k in 1:dim.m)
  {
    for(j in 1:n.comp)
    {
      result[k, j] = ind.limit.y.pos(comp.dist[k,j], params.list[[k]][[j]], limit = limit[k])
    }
  }

  return(result)
}

#' A function to calculate the limited mean of a matrix of given distributions (with zero inflation), by component and by dimension.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param limit A vector of limit values by dimension.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of limited mean values by dimension and by component.
#'
#' @keywords internal
#'
# #' @export dim.comp.limit.y
dim.comp.limit.y = function(comp.dist, zero.prob, params.list, limit)
{
  temp = dim.comp.limit.y.pos(comp.dist, params.list, limit)
  result = (1-zero.prob)*temp

  return(result)
}

#' Predict the limited mean of y, given a fixed covariate vector X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param limit A vector of limit values by dimension.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A vector of limited mean values by dimension.
#'
#' @keywords internal
#'
# #' @export ind.predict.limit
ind.predict.limit = function(X, alpha, comp.dist, zero.prob, params.list, limit)
{
  n.comp = nrow(alpha)
  dim.m = nrow(comp.dist)

  weighting = exp(gate.logit(X, alpha))
  temp = dim.comp.limit.y(comp.dist, zero.prob, params.list, limit)
  result = temp%*%t(weighting)

  return(t(result))
}

#' Predict the limited mean of y, given a fixed covariate matrix X, a model and a vector of limits by dimension.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param limit A vector of limit to apply to each dimension of y.
#'
#' @seealso \code{\link{LRMoE.fit}}, \code{\link[actuar]{GammaSupp}}.
#'
#' @return A matrix of limited mean values by observation and by dimension.
#'         Calculation is done for severity distributions only. \code{NA} values are returned for frequency distributions.
#'
#' @rawNamespace S3method(predict, limit)
#'
#' @export predict.limit
predict.limit = function(X, alpha, comp.dist, zero.prob, params.list, limit)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)

  if(is.null(X.size))
  {
    return( ind.predict.limit(X, alpha, comp.dist, zero.prob, params.list, limit)  )
  }

  result = array(0, dim = c(X.size, dim.m))
  # for(i in 1:(X.size) )
  # {
  #   result[i,] = ind.predict.limit(X[i,], alpha, comp.dist, zero.prob, params.list, limit)
  # }
  result = apply(X, MARGIN = 1, FUN = ind.predict.limit,
                 alpha = alpha, comp.dist = comp.dist, params.list = params.list, zero.prob = zero.prob, limit = limit)

  return(t(t(result)))
}

#' Predict the excess mean of y, given a fixed covariate matrix X, a model and a vector of limits by dimension.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param limit A vector of limit to apply to each dimension of y.
#'
#' @seealso \code{\link{LRMoE.fit}}, \code{\link{predict.limit}}.
#'
#' @return A matrix of excess mean values by observation and by dimension.
#'         This is equan to \code{predict.mean - predict.limit}.
#'
#' @rawNamespace S3method(predict, excess)
#'
#' @export predict.excess
predict.excess = function(X, alpha, comp.dist, zero.prob, params.list, limit)
{
  result = predict.mean(X, alpha, comp.dist, zero.prob, params.list) - predict.limit(X, alpha, comp.dist, zero.prob, params.list, limit)

  return(result)
}





