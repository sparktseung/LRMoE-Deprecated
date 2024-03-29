## Expert Function: Weibull
#' Expert Function: Weibull.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param shape.k A vector of length \code{g}: Weibull shape parameters.
#' @param scale.lambda A vector of length \code{g}: Weibull scale parameters.
#' @return A list of matrices of expert loglikelihood for Weibull.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[stats]{Weibull}}.
#'
#' @importFrom stats pweibull dweibull
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
# #' @export expert.weibull
expert.weibull = function(tl, yl, yu, tu, g = 1, shape.k, scale.lambda)
{
  # Initialization: return value are N * g matrices
  expert.weibull.ll=expert.weibull.tn=expert.weibull.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=stats::pweibull(yu[censor.idx],shape = shape.k[j],scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=stats::pweibull(yl[censor.idx],shape = shape.k[j],scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.weibull.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.weibull.ll[!censor.idx,j]=stats::dweibull(yu[!censor.idx],shape = shape.k[j],scale = scale.lambda[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=stats::pweibull(tu,shape = shape.k[j],scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=stats::pweibull(tl,shape = shape.k[j],scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.weibull.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.weibull.tn[no.trunc.idx,j] = stats::dweibull(tu[no.trunc.idx],shape = shape.k[j],scale = scale.lambda[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.weibull.ll[zero.idx,j]=(-Inf)
    expert.weibull.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.weibull.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.weibull.tn[!no.trunc.idx,j])
    expert.weibull.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.weibull.ll=expert.weibull.ll, expert.weibull.tn=expert.weibull.tn, expert.weibull.tn.bar=expert.weibull.tn.bar)
}
