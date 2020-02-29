## Expert Function: Gamma Count
#' Expert Function: Gamma Count.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param m A vector of length \code{g}: Gamma Count shape parameters.
#' @param s A vector of length \code{g}: Gamma Count dispersion parameters.
#' @return A list of matrices of expert loglikelihood for Gamma Count.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[rmutil]{GammaCount}}.
#'
#' @importFrom rmutil pgammacount dgammacount
#' @importFrom copula log1mexp
#'
#' @keywords internal
#'
#' @export expert.gammacount
expert.gammacount = function(tl, yl, yu, tu, g = 1, m, s)
{
  # Initialization: return value are N * g matrices
  expert.gammacount.ll=expert.gammacount.tn=expert.gammacount.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)

    # prob.log.yu = log(rmutil::pgammacount(yu, m[j], s[j]))
    prob.log.yu = pgammacount.new(yu[censor.idx], m = m[j], s = s[j], log = TRUE)
    # prob.log.yl = log(rmutil::pgammacount(yl, m[j], s[j]))
    # prob.log.yl = pgammacount.new(yl[censor.idx], m = m[j], s = s[j], log = TRUE)
    prob.log.yl = pgammacount.new(ceiling(yl[censor.idx])-1, m = m[j], s = s[j], log = TRUE)

    # Compute loglikelihood for expert j, first for y
    expert.gammacount.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.gammacount.ll[!censor.idx,j]= dgammacount.new(yl[!censor.idx], m = m[j], s = s[j], log = TRUE)

    # Compute loglikelihood for expert j, then for truncation limits t
    # prob.log.tu = log(rmutil::pgammacount(tu, m[j], s[j]))
    prob.log.tu = pgammacount.new(tu, m = m[j], s = s[j], log = TRUE)
    # prob.log.tl = log(rmutil::pgammacount(tl, m[j], s[j]))
    # prob.log.tl = ifelse(tl==0, -Inf, log(rmutil::pgammacount(tl, m[j], s[j])))
    # prob.log.tl = pgammacount.new(tl, m = m[j], s = s[j], log = TRUE)
    prob.log.tl = pgammacount.new(ceiling(tl)-1, m = m[j], s = s[j], log = TRUE)

    # Normalizing factor for truncation limits, in log
    expert.gammacount.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    expert.gammacount.tn[no.trunc.idx,j] = dgammacount.new(tu[no.trunc.idx], m =  m[j], s = s[j], log = FALSE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    # zero.idx = (tu==0)
    # expert.gammacount.ll[zero.idx,j]=(-Inf)
    # expert.gammacount.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.gammacount.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.gammacount.tn[!no.trunc.idx,j])
    expert.gammacount.tn.bar[no.trunc.idx,j] = log1mexp(-expert.gammacount.tn[!no.trunc.idx,j])
  }
  # Return values
  list(expert.gammacount.ll=expert.gammacount.ll, expert.gammacount.tn=expert.gammacount.tn, expert.gammacount.tn.bar=expert.gammacount.tn.bar)
}
