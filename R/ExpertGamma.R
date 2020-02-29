## Expert Function: Gamma
#' Expert Function: Gamma.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param m A vector of length \code{g}: Gamma shape parameters.
#' @param theta A vector of length \code{g}: Gamma scale parameters.
#' @return A list of matrices of expert loglikelihood for Gamma.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[stats]{GammaDist}}.
#'
#' @importFrom stats pgamma dgamma
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export expert.gamma
expert.gamma = function(tl, yl, yu, tu, g = 1, m, theta)
{
  # Initialization: return value are N * g matrices
  expert.gamma.ll=expert.gamma.tn=expert.gamma.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=pgamma(yu[censor.idx],shape=m[j],scale=theta[j],log.p=TRUE)
    prob.log.yl=pgamma(yl[censor.idx],shape=m[j],scale=theta[j],log.p=TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.gamma.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.gamma.ll[!censor.idx,j]=dgamma(yu[!censor.idx],shape=m[j],scale=theta[j],log=TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=pgamma(tu,shape=m[j],scale=theta[j],log.p=TRUE)
    prob.log.tl=pgamma(tl,shape=m[j],scale=theta[j],log.p=TRUE)

    # Normalizing factor for truncation limits, in log
    expert.gamma.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.gamma.tn[no.trunc.idx,j] = dgamma(tu[no.trunc.idx],shape=m[j],scale=theta[j],log=TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.gamma.ll[zero.idx,j]=(-Inf)
    expert.gamma.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.gamma.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.gamma.tn[!no.trunc.idx,j])
    expert.gamma.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.gamma.ll=expert.gamma.ll, expert.gamma.tn=expert.gamma.tn, expert.gamma.tn.bar=expert.gamma.tn.bar)
}
