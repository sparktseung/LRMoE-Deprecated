## Expert Function: Poisson
#' Expert Function: Poisson.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param mean.theta A vector of length \code{g}: Poisson mean parameters.
#' @return A list of matrices of expert loglikelihood for Poisson.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[stats]{Poisson}}.
#'
#' @importFrom stats ppois dpois
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export expert.poisson
expert.poisson = function(tl, yl, yu, tu, g = 1, mean.theta)
{
  # Initialization: return value are N * g matrices
  expert.poisson.ll=expert.poisson.tn=expert.poisson.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=ppois(yu[censor.idx],mean.theta[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=ppois(ceiling(yl[censor.idx])-1,mean.theta[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.poisson.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.poisson.ll[!censor.idx,j]=dpois(yu[!censor.idx],mean.theta[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=ppois(tu,mean.theta[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=ppois(ceiling(tl)-1,mean.theta[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.poisson.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.poisson.tn[no.trunc.idx,j] = dpois(tu[no.trunc.idx],mean.theta[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case: The following code treating zeros is NOT applicable for frequency distributions!!!
    # zero.idx = (tu==0)
    # expert.poisson.ll[zero.idx,j]=(-Inf)
    # expert.poisson.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.poisson.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.poisson.tn[!no.trunc.idx,j])
    expert.poisson.tn.bar[no.trunc.idx,j] = log1mexp(-expert.poisson.tn[no.trunc.idx,j])
    # expert.poisson.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.poisson.ll=expert.poisson.ll, expert.poisson.tn=expert.poisson.tn, expert.poisson.tn.bar=expert.poisson.tn.bar)
}
