#' Count the number of logit regression coefficients.
#'
#' @param alpha A g*P matrix of logit regression coefficients.
#'
#' @return The number of logit regression coefficients.
#'
#'
#' @keywords internal
#'
#' @export count.alpha
count.alpha = function(alpha)
{
  return( ncol(alpha)*(nrow(alpha)-1) )
}



#' Count the number of zero probability masses.
#'
#' @param comp.dist A d*g matrix of strings describing the component distributions.
#'
#' @return The number of zero probability masses.
#'
#'
#' @keywords internal
#'
#' @export count.zero
count.zero = function(comp.dist)
{
  result = 0
  for(k in 1:(nrow(comp.dist)))
  {
    for(j in 1:(ncol(comp.dist)))
    {
      if(substr(comp.dist[k,j], 1, 3)=="ZI-")
      {
        result = result + 1
      }
    }
  }
  return(result)
}



#' Count the number of parameters in the positive part.
#'
#' @param comp.dist A d*g matrix of strings describing the component distributions.
#'
#' @return The number of parameters in the positive part.
#'
#'
#' @keywords internal
#'
#' @export count.pos.params
count.pos.params = function(comp.dist)
{
  result = 0

  for(k in 1:(nrow(comp.dist)))
  {
    for(j in 1:(ncol(comp.dist)))
    {
      switch (comp.dist[k,j],
              # Severity distributions & their zero-inflation
              "gamma"       = {result = result + 2},
              "ZI-gamma"    = {result = result + 2},
              "invgauss"    = {result = result + 2},
              "ZI-invgauss" = {result = result + 2},
              "lnorm"       = {result = result + 2},
              "ZI-lnorm"    = {result = result + 2},
              "weibull"     = {result = result + 2},
              "ZI-weibull"  = {result = result + 2},
              "burr"        = {result = result + 3},
              "ZI-burr"     = {result = result + 3},
              # Frequency distributions & their zero-inflation
              "poisson"     = {result = result + 1},
              "ZI-poisson"  = {result = result + 1},
              "nbinom"      = {result = result + 2},
              "ZI-nbinom"   = {result = result + 2},
              "gammacount"  = {result = result + 2},
              "ZI-gammacount"  = {result = result + 2}
      )
    }
  }
  return(result)
}
