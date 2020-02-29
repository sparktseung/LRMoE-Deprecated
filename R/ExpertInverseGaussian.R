## Expert Function: Inverse Gaussian
#' Expert Function: Inverse Gaussian.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param mean.mu A vector of length \code{g}: Inverse Gaussian mean parameters.
#' @param shape.lambda A vector of length \code{g}: Inverse Gaussian shape parameters.
#' @return A list of matrices of expert loglikelihood for Inverse Gaussian.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[statmod]{invgauss}}.
#'
#' @importFrom statmod pinvgauss dinvgauss
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export expert.invgauss
expert.invgauss = function(tl, yl, yu, tu, g = 1, mean.mu, shape.lamba)
{
  # Initialization: return value are N * g matrices
  expert.invgauss.ll=expert.invgauss.tn=expert.invgauss.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=statmod::pinvgauss(yu[censor.idx],mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=statmod::pinvgauss(yl[censor.idx],mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.invgauss.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.invgauss.ll[!censor.idx,j]=statmod::dinvgauss(yu[!censor.idx],mean = mean.mu[j],shape = shape.lamba[j], log = TRUE) # exact likelihood

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=statmod::pinvgauss(tu,mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=statmod::pinvgauss(tl,mean = mean.mu[j],shape = shape.lamba[j], lower.tail = TRUE, log.p = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.invgauss.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.invgauss.tn[no.trunc.idx,j] = statmod::dinvgauss(tu[no.trunc.idx],mean = mean.mu[j],shape = shape.lamba[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.invgauss.ll[zero.idx,j]=(-Inf)
    expert.invgauss.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.invgauss.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.invgauss.tn[!no.trunc.idx,j])
    expert.invgauss.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.invgauss.ll=expert.invgauss.ll, expert.invgauss.tn=expert.invgauss.tn, expert.invgauss.tn.bar=expert.invgauss.tn.bar)
}
