## Positive part of Expert functions
#' A switch to calculate the loglikelihood of expert functions.
#'
#' @param ind.dist A string which indicates the expert function.
#' \itemize{
#'     \item \code{gamma}: Gamma
#'     \item \code{lnorm}: Log Normal
#'     \item \code{invgauss}: Inverse Gaussian
#'     \item \code{weibull}: Weibull
#'     \item \code{burr}: Burr
#'     \item \code{poisson}: Poisson
#'     \item \code{nbinom}: Negative Binomial
#'     \item \code{binom}: Binomial
#'     \item \code{gammacount}: Gamma Count
#'     \item \code{ZI-root}: Zero-inflated versions of the distributions above, e.g. \code{ZI-gamma}.
#' }
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param params A vector of parameters for \code{ind.dist}.
#'
#' @return A list of matrices of expert loglikelihood for the chosen distribution \code{ind.dist}.
#' \itemize{
#'     \item \code{expert.ind.dist.ll}: An N*1 matrix of loglikelihood: density for exact observation; probability for censored observation.
#'     \item \code{expert.ind.dist.tn}: An N*1 matrix of loglikelihood within truncation limits.
#'     \item \code{expert.ind.dist.tn.bar}: An N*1 matrix of loglikelihood outside truncation limits, which is \code{log(1-exp(expert.ind.dist.tn))}.
#' }
#'
#' @seealso \code{\link{expert.gamma}}, \code{\link{expert.lognormal}}, \code{\link{expert.invgauss}}, \code{\link{expert.weibull}}, \code{\link{expert.burr}},
#'          \code{\link{expert.poisson}}, \code{\link{expert.nbinom}}, \code{\link{expert.gammacount}}
#'
#' @keywords internal
#'
# #' @export pos.expert.loglik.calc
pos.expert.loglik.calc = function(ind.dist, tl, yl, yu, tu, params)
{
  temp = NULL
  switch (ind.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = expert.gamma(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-gamma"    = {temp = expert.gamma(tl, yl, yu, tu, 1, params[1], params[2])},
          "invgauss"    = {temp = expert.invgauss(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-invgauss" = {temp = expert.invgauss(tl, yl, yu, tu, 1, params[1], params[2])},
          "lnorm"       = {temp = expert.lognormal(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-lnorm"    = {temp = expert.lognormal(tl, yl, yu, tu, 1, params[1], params[2])},
          "weibull"     = {temp = expert.weibull(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-weibull"  = {temp = expert.weibull(tl, yl, yu, tu, 1, params[1], params[2])},
          "burr"        = {temp = expert.burr(tl, yl, yu, tu, 1, params[1], params[2], params[3])},
          "ZI-burr"     = {temp = expert.burr(tl, yl, yu, tu, 1, params[1], params[2], params[3])},
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = expert.poisson(tl, yl, yu, tu, 1, params[1])},
          "ZI-poisson"  = {temp = expert.poisson(tl, yl, yu, tu, 1, params[1])},
          "nbinom"      = {temp = expert.nbinom(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-nbinom"   = {temp = expert.nbinom(tl, yl, yu, tu, 1, params[1], params[2])},
          "binom"       = {temp = expert.binom(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-binom"    = {temp = expert.binom(tl, yl, yu, tu, 1, params[1], params[2])},
          "gammacount"  = {temp = expert.gammacount(tl, yl, yu, tu, 1, params[1], params[2])},
          "ZI-gammacount"  = {temp = expert.gammacount(tl, yl, yu, tu, 1, params[1], params[2])},
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
