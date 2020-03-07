## Expert Function: Nebative Binomial
#' Expert Function: Negative Binomial.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param size.n A vector of length \code{g}: Negative Binomial size.n parameters.
#' @param prob.p A vector of length \code{g}: Negative Binomial prob.p parameters.
#' @return A list of matrices of expert loglikelihood for Negative Binomial.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[stats]{NegBinomial}}.
#'
#' @importFrom stats pnbinom dnbinom
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
# #' @export expert.nbinom
expert.nbinom = function(tl, yl, yu, tu, g = 1, size.n, prob.p)
{
  # Initialization: return value are N * g matrices
  expert.nbinom.ll=expert.nbinom.tn=expert.nbinom.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=pnbinom(yu[censor.idx], size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=pnbinom(ceiling(yl[censor.idx])-1, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.nbinom.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.nbinom.ll[!censor.idx,j]=dnbinom(yu[!censor.idx], size = size.n[j], prob = prob.p[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=pnbinom(tu, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=pnbinom(ceiling(tl)-1, size = size.n[j], prob = prob.p[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.nbinom.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.nbinom.tn[no.trunc.idx,j] = dnbinom(tu[no.trunc.idx], size = size.n[j], prob = prob.p[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case: The following code treating zeros is NOT applicable for frequency distributions!!!
    # zero.idx = (tu==0)
    # expert.nbinom.ll[zero.idx,j]=(-Inf)
    # expert.nbinom.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.nbinom.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.nbinom.tn[!no.trunc.idx,j])
    expert.nbinom.tn.bar[no.trunc.idx,j] = log1mexp(-expert.nbinom.tn[no.trunc.idx,j])
    # expert.nbinom.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.nbinom.ll=expert.nbinom.ll, expert.nbinom.tn=expert.nbinom.tn, expert.nbinom.tn.bar=expert.nbinom.tn.bar)
}
