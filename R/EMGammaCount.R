## ECM algorithm of Gamma Count expert
#' ECM: M-Step for Gamma Count expert.
#'
#' @seealso \code{\link{comp.kj.params.m.recur}}
#'
#' @importFrom stats dnbinom
#'
#' @keywords internal
#'
#' @export gammacount.params.m.recur
gammacount.params.m.recur = function(gammacount.params.old,
                                     tl, yl, yu, tu,
                                     expert.ll, expert.tn, expert.tn.bar,
                                     z.e.obs, z.e.lat, k.e,
                                     penalty, hyper.params)
{
  # Constants
  sample.size.n = length(tl)
  # Take old parameters
  shape.m = gammacount.params.old[1]
  disp.s = gammacount.params.old[2]
  # Take hyper parameters
  hyper.m.1 = hyper.params[1]
  hyper.m.2 = hyper.params[2]
  hyper.s.1 = hyper.params[3]
  hyper.s.2 = hyper.params[4]
  # Value to return
  gammacount.params.new = gammacount.params.old

  # Use brute-force optimization, due to complex form of pdf/cdf of gammacount
  Q.params = function(params.new,
                      tl, yl, yu, tu,
                      z.e.obs, z.e.lat, k.e,
                      penalty, hyper.m.1, hyper.m.2, hyper.s.1, hyper.s.2)
  {
    censor.idx = (yl!=yu)
    sample.size.n = length(tl)

    expert.gammacount.ll=expert.gammacount.tn=expert.gammacount.tn.bar=array(-Inf, dim=c(sample.size.n,1))

    shape.m.new = params.new[1]
    disp.s.new = params.new[2]

    prob.log.yu = ifelse(yu[censor.idx]==Inf, 0, pgammacount.new(yu[censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE))
    prob.log.yl = ifelse(yl[censor.idx]==0, -Inf, pgammacount.new(ceiling(yl[censor.idx])-1, m = shape.m.new, s = disp.s.new, log = TRUE))

    expert.gammacount.ll[censor.idx,1]=prob.log.yu+log1mexp(prob.log.yu-prob.log.yl) # likelihood of censored interval: some easy algebra
    expert.gammacount.ll[!censor.idx,1]= dgammacount.new(yl[!censor.idx], m = shape.m.new, s = disp.s.new, log = TRUE)

    prob.log.tu = ifelse(tu==Inf, 0, pgammacount.new(tu, m = shape.m.new, s = disp.s.new, log = TRUE))
    prob.log.tl = ifelse(tl==0, -Inf, pgammacount.new(ceiling(tl)-1, m = shape.m.new, s = disp.s.new, log = TRUE))

    expert.gammacount.tn[,1]=prob.log.tu+log1mexp(prob.log.tu-prob.log.tl)
    # Deal with no truncation case
    no.trunc.idx = (tl==tu)

    # expert.gammacount.tn[no.trunc.idx,1] = rmutil::dgammacount(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)
    expert.gammacount.tn[no.trunc.idx,1] = dgammacount.new(tu[no.trunc.idx], m = shape.m.new, s = disp.s.new, log = TRUE)

    expert.gammacount.tn.bar[!no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[!no.trunc.idx,1])
    expert.gammacount.tn.bar[no.trunc.idx,1] = log1mexp(-expert.gammacount.tn[no.trunc.idx,1])

    # result = apply(sweep(matrix(z.e.obs), 1, matrix(expert.gammacount.ll), FUN = "*", check.margin = FALSE), 2, FUN = "sum") +
    #   apply(sweep(matrix(k.e), 2,
    #               sweep(matrix(z.e.lat), 2, matrix(expert.gammacount.tn), FUN = "*", check.margin = FALSE),
    #               FUN = "*", check.margin = FALSE),
    #         2, FUN = "sum")
    result = sum(z.e.obs*expert.gammacount.ll) + sum(k.e*z.e.lat*expert.gammacount.tn)
    if(penalty==TRUE){
      result = result + (hyper.m.1-1)*log(shape.m.new) - shape.m.new/hyper.m.2 + (hyper.s.1-1)*log(disp.s.new) - disp.s.new/hyper.s.2
    }

    return(result * (-1))
  }

  # pos.idx = (yu!=0)
  pos.idx = rep(TRUE, length(yu))

  temp.params = optim(par = gammacount.params.old, fn = Q.params,
                      tl = tl[pos.idx], yl = yl[pos.idx], yu = yu[pos.idx], tu = tu[pos.idx],
                      z.e.obs = z.e.obs[pos.idx], z.e.lat = z.e.lat[pos.idx], k.e = k.e[pos.idx],
                      penalty = penalty, hyper.m.1 = hyper.m.1, hyper.m.2 = hyper.m.2, hyper.s.1 = hyper.s.1, hyper.s.2 = hyper.s.2,
                      method = "L-BFGS-B", lower = 0.5*gammacount.params.old, upper = 2*gammacount.params.old)$par

  gammacount.params.new = temp.params


  return(gammacount.params.new)
}
