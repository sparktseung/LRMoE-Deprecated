## Expert Function: Burr
#' Expert Function: Burr.
#'
#' @param tl A vector of length N: lower bounds of truncation.
#' @param yl A vector of length N: lower bounds of censoring.
#' @param yu A vector of length N: upper bounds of censoring.
#' @param tu A vector of length N: upper bounds of truncation.
#' @param g A numeric which is set to 1. See Note.
#' @param shape1.k A vector of length \code{g}: Burr shape1 parameters.
#' @param shape2.c A vector of length \code{g}: Burr shape2 parameters.
#' @param scale.lambda A vector of length \code{g}: Burr scale parameters.
#' @return A list of matrices of expert loglikelihood for Burr.
#'
#' @section Note:
#' \code{g} is a legacy input. It is always set to 1, because this expert function is called within function \code{pos.expert.loglik.calc}.
#'
#' @seealso \code{\link{pos.expert.loglik.calc}}, \code{\link[actuar]{Burr}}.
#'
#' @importFrom actuar pburr dburr
#' @importFrom copula log1mexp log1pexp
#'
#' @keywords internal
#'
#' @export expert.burr
expert.burr = function(tl, yl, yu, tu, g = 1, shape1.k, shape2.c, scale.lambda)
{
  # Initialization: return value are N * g matrices
  expert.burr.ll=expert.burr.tn=expert.burr.tn.bar=array(-Inf, dim=c(length(yu),g))

  # for (j in 1:g) # Loop trough each expert component
  j = 1
  {
    # Find indexes of unequal yl & yu: Not exact observations, but censored
    censor.idx=(yl!=yu)
    prob.log.yu=actuar::pburr(yu[censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.yl=actuar::pburr(yl[censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    # prob.log.yu = log1mexp( shape1.k[j] * log(1+(yu[censor.idx]/scale.lambda[j])^(shape2.c[j])) )
    # prob.log.yl = log1mexp( shape1.k[j] * log(1+(yl[censor.idx]/scale.lambda[j])^(shape2.c[j])) )

    # Compute loglikelihood for expert j, first for y
    expert.burr.ll[censor.idx,j]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra

    ###################################################################
    # Deal with numerical underflow: prob.log.yu and prob.log.yl can both be -Inf
    NA.idx = which(is.na(expert.burr.ll[,j]))
    expert.burr.ll[NA.idx, j] = -Inf
    ###################################################################

    expert.burr.ll[!censor.idx,j]=actuar::dburr(yu[!censor.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], log = TRUE) # exact likelihood
    # expert.burr.ll[!censor.idx,j] = log(shape1.k[j]*shape2.c[j]) + (shape2.c[j]-1) * log(yu[!censor.idx]) - shape2.c[j] * log(scale.lambda[j]) +
    #                                 (-shape1.k[j]-1) * log1pexp(shape2.c[j] * log(yu[!censor.idx]/scale.lambda[j]) )

    # Compute loglikelihood for expert j, then for truncation limits t
    prob.log.tu=actuar::pburr(tu,shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    prob.log.tl=actuar::pburr(tl,shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], lower.tail = TRUE, log.p = TRUE)
    # prob.log.tu = log1mexp( shape1.k[j] * log(1+(tu/scale.lambda[j])^(shape2.c[j])) )
    # prob.log.tl = log1mexp( shape1.k[j] * log(1+(tl/scale.lambda[j])^(shape2.c[j])) )

    # Normalizing factor for truncation limits, in log
    expert.burr.tn[,j]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)

    ###################################################################
    # Deal with numerical underflow: prob.log.tu and prob.log.tl can both be -Inf
    NA.idx = which(is.na(expert.burr.tn[,j]))
    expert.burr.tn[NA.idx, j] = -Inf
    ###################################################################

    ###################################################################
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)
    # expert.burr.tn[no.trunc.idx,j] = log(shape1.k[j]*shape2.c[j]) + (shape2.c[j]-1) * log(tu[no.trunc.idx]) - shape2.c[j] * log(scale.lambda[j]) +
    #   (-shape1.k[j]-1) * log1pexp(shape2.c[j] * log(tu[no.trunc.idx]/scale.lambda[j]) )
    expert.burr.tn[no.trunc.idx,j] = actuar::dburr(tu[no.trunc.idx],shape1 = shape1.k[j],shape2 = shape2.c[j], scale = scale.lambda[j], log = TRUE)
    ###################################################################

    ###################################################################
    # Deal with exact zero case
    zero.idx = (tu==0)
    expert.burr.ll[zero.idx,j]=(-Inf)
    expert.burr.tn[zero.idx,j]=(-Inf)
    ###################################################################

    # Log of Pr(outside of truncation interval)
    expert.burr.tn.bar[!no.trunc.idx,j] = log1mexp(-expert.burr.tn[!no.trunc.idx,j])
    expert.burr.tn.bar[no.trunc.idx,j] = 0
  }
  # Return values
  list(expert.burr.ll=expert.burr.ll, expert.burr.tn=expert.burr.tn, expert.burr.tn.bar=expert.burr.tn.bar)
}
