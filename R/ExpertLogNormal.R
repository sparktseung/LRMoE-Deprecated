## Expert Function: Log Normal
#' Expert Function: Log Normal.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param meanlog A vector of length \code{g}: Lognormal mean parameters.
#' @param sdlog A vector of length \code{g}: Lognormal sd parameters.
#' @return A list of matrices of expert loglikelihood for Log Normal.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[stats]{Lognormal}}.
#'
#' @importFrom stats plnorm dlnorm
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
# #' @export expert.lognormal
expert.lognormal = function(tl, yl, yu, tu, g = 1, meanlog, sdlog)
{
  # Initialization: return value are N * g matrices
  expert.lognormal.ll=expert.lognormal.tn=expert.lognormal.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=plnorm(yu[censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)
    prob.log.yl=plnorm(yl[censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.lognormal.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.lognormal.ll[!censor.idx,j]=dlnorm(yu[!censor.idx],meanlog = meanlog[j],sdlog = sdlog[j],log=TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=plnorm(tu,meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)
    prob.log.tl=plnorm(tl,meanlog = meanlog[j],sdlog = sdlog[j],log.p=TRUE)

    # Normalizing factor for truncation limits, in log
    expert.lognormal.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.lognormal.tn[no.trunc.idx,j] = dlnorm(tu[no.trunc.idx],meanlog = meanlog[j],sdlog = sdlog[j],log=TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.lognormal.ll[zero.idx,j]=(-Inf)
    expert.lognormal.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.lognormal.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.lognormal.tn[!no.trunc.idx,j])
    expert.lognormal.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.lognormal.ll=expert.lognormal.ll, expert.lognormal.tn=expert.lognormal.tn, expert.lognormal.tn.bar=expert.lognormal.tn.bar)
}
