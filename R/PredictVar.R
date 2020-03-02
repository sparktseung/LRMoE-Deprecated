#' A switch to calculate the variance of a matrix of given distributions (positive part only), by dimension and by component.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of variance by dimension and by component.
#'
#' @keywords internal
#'
#' @export dim.comp.var.y.pos
dim.comp.var.y.pos = function(comp.dist, params.list)
{
  n.comp = ncol(comp.dist)
  dim.m = nrow(comp.dist)

  result = matrix(0, nrow = dim.m, ncol = n.comp)

  for(k in 1:dim.m)
  {
    for(j in 1:n.comp)
    {
      result[k, j] = ind.var.y.pos(comp.dist[k,j], params.list[[k]][[j]])
    }
  }

  return(result)
}


#' A function to calculate the variance of a matrix of given distributions (with zero inflation), by dimension and by component.
#'
#' @param comp.dist A d*g matrix of strings, describing component distributions by dimension and by component.
#' @param zero.prob A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#' @param paramas.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of variance by dimension and by component.
#'
#' @keywords internal
#'
#' @export dim.comp.var.y
dim.comp.var.y = function(comp.dist, zero.prob, params.list)
{
  cond.mean = dim.comp.mean.y.pos(comp.dist, params.list)
  cond.var = dim.comp.var.y.pos(comp.dist, params.list)

  result = (1-zero.prob)*cond.var + zero.prob*(1-zero.prob)*cond.mean^2

  return(result)
}


#' Predict the variance of y, given a fixed covariate vector X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A vector of variance by dimension.
#'
#' @keywords internal
#'
#' @export ind.predict.var
ind.predict.var = function(X, alpha, comp.dist, zero.prob, params.list)
{
  n.comp = nrow(alpha)
  dim.m = nrow(comp.dist)

  weighting = exp(gate.logit(X, alpha))

  # Var of conditional expectation
  cond.mean = dim.comp.mean.y(comp.dist, zero.prob, params.list)
  cond.mean.squared = cond.mean^2
  grand.mean = cond.mean%*%t(weighting)
  var.cond.mean = cond.mean.squared%*%t(weighting) - (grand.mean)^2

  # Mean of conditional variance
  cond.var = dim.comp.var.y(comp.dist, zero.prob, params.list)
  ave.cond.var = cond.var%*%t(weighting)

  return(t(var.cond.mean + ave.cond.var))
}


#' Predict the variance of y, given a fixed covariate matrix X and a model.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A matrix of variance by observation and by dimension.
#'
#' @export predict.var
predict.var = function(X, alpha, comp.dist, zero.prob, params.list)
{
  X.size = nrow(X)
  dim.m = nrow(comp.dist)

  if(is.null(X.size))
  {
    return( ind.predict.var(X, alpha, comp.dist, zero.prob, params.list)  )
  }

  result = array(0, dim = c(X.size, dim.m))
  # for(i in 1:(X.size) )
  # {
  #   result[i,] = ind.predict.var(X[i,], alpha, comp.dist, zero.prob, params.list)
  # }

  result = apply(X, MARGIN = 1, FUN = ind.predict.var,
                 alpha = alpha, comp.dist = comp.dist, zero.prob = zero.prob, params.list = params.list)

  return(t(result))
}




