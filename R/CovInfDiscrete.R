#' Investigate the influence of X (discrete covariate) on Y.
#'
#' @param X A matrix of covariates.
#' @param idx Column index of the discrete covariate of interest. There is no need to include the intercept.
#' @param response A string vector indicating what metrics of Y to calculate. Currently supported calculation include:
#' \itemize{
#'    \item "Mean" : Mean of response
#'    \item "SD" : Variance of response
#'    \item "VAR990", "CTE990" : Value-at-Risk and Conditional Tail Expectation of response.
#'           The last three digits indicate the desired probability level, e.g. "990" = 99.0%.
#' }
#' @param dim The dimension of the response under consideration.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @seealso \code{\link{LRMoE.fit}}.
#'
#' @return A data frame containing all results.
#'
#' @export covinf.discrete
covinf.discrete = function(X, idx, response = c("Mean", "VAR990", "CTE990"), dim = 1,
                           alpha, comp.dist, zero.prob, params.list)
{
  # Save a copy of X
  X.temp = X
  # Result to return
  df = data.frame(matrix(ncol=length(response),nrow=1+length(idx),
                         dimnames=list(colnames(X[,c(1,idx)]), response)))

  # Base case
  X.temp[,idx] = 0
  for(j in 1:length(response)){
    df[1,j] = mean(
      switch (substr(response[j], 1, 2),
              "Me" = LRMoE::predict.mean(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
              "SD" = LRMoE::predict.var(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
              "VA" = LRMoE::predict.quantile(X.temp, alpha, comp.dist, zero.prob, params.list,
                                             rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
              "CT" = LRMoE::predict.cte(X.temp, alpha, comp.dist, zero.prob, params.list,
                                        rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
              # Error
              stop("Invalid input!")
      )
    )
  }

  # Other cases
  for(k in 1:length(idx))
  {
    X.temp[,idx[k]] = 1
    if(k>1){X.temp[,idx[k-1]] = 0}
    for(j in 1:length(response)){
      df[k+1,j] = mean(
        switch (substr(response[j], 1, 2),
                "Me" = LRMoE::predict.mean(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "SD" = LRMoE::predict.var(X.temp, alpha, comp.dist, zero.prob, params.list)[,dim],
                "VA" = LRMoE::predict.quantile(X.temp, alpha, comp.dist, zero.prob, params.list,
                                               rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                "CT" = LRMoE::predict.cte(X.temp, alpha, comp.dist, zero.prob, params.list,
                                          rep(strtoi(substr(response[j], 4, 6))/1000, nrow(comp.dist)) )[,dim],
                # Error
                stop("Invalid input!")
        )
      )
    }

  }

  return(df)

}
