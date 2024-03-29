## ECM: M-Step for updating parameters of positive part of component experts.
#' ECM: M-Step for updating parameters of positive part of component experts.
#'
#' @param comp.kj.dist A string which indicates component distribution. See also \code{\link{pos.expert.loglik.calc}}.
#' @param comp.kj.params.old A vector of numerics: old parameter values for \code{comp.kj.dist}.
#' @param tl.k,yl.k,yu.k,tu.k Numerical vectors of length N. See also \code{\link{pos.expert.loglik.calc}}.
#' @param comp.kj.pos.expert.ll,comp.kj.pos.expert.tn,comp.kj.pos.expert.tn.bar See \code{\link{pos.expert.loglik.calc}}.
#' @param z.e.obs,z.e.lat,k.e A numerical vector returned by \code{\link{comp.zkz.e.recur}}.
#' @param penalty,hyper.params.kj See \code{\link{expert.loglik.pen.dim.comp}}.
#'
#' @return Updated parameter values.
#'
#' @keywords internal
#'
# #' @export comp.kj.params.m.recur
comp.kj.params.m.recur = function(comp.kj.dist, comp.kj.params.old,
                                  tl.k, yl.k, yu.k, tu.k,
                                  comp.kj.pos.expert.ll, comp.kj.pos.expert.tn, comp.kj.pos.expert.tn.bar,
                                  z.e.obs, z.e.lat, k.e,
                                  penalty, hyper.params.kj)
{
  temp = NULL
  switch (comp.kj.dist,
          # Severity distributions & their zero-inflation
          "gamma"       = {temp = gamma.params.m.recur(gamma.params.old = comp.kj.params.old,
                                                       tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                       expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                       z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                       penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-gamma"    = {temp = gamma.params.m.recur(gamma.params.old = comp.kj.params.old,
                                                       tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                       expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                       z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                       penalty = penalty, hyper.params = hyper.params.kj) },
          "invgauss"    = {temp = invgauss.params.m.recur(invgauss.params.old = comp.kj.params.old,
                                                          tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                          expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                          z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                          penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-invgauss" = {temp = invgauss.params.m.recur(invgauss.params.old = comp.kj.params.old,
                                                          tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                          expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                          z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                          penalty = penalty, hyper.params = hyper.params.kj) },
          "lnorm"       = {temp = lnorm.params.m.recur(lnorm.params.old = comp.kj.params.old,
                                                       tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                       expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                       z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                       penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-lnorm"    = {temp = lnorm.params.m.recur(lnorm.params.old = comp.kj.params.old,
                                                       tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                       expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                       z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                       penalty = penalty, hyper.params = hyper.params.kj) },
          "weibull"     = {temp = weibull.params.m.recur(weibull.params.old = comp.kj.params.old,
                                                         tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                         expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                         z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                         penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-weibull"  = {temp = weibull.params.m.recur(weibull.params.old = comp.kj.params.old,
                                                         tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                         expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                         z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                         penalty = penalty, hyper.params = hyper.params.kj) },
          "burr"        = {temp = burr.params.m.recur(burr.params.old = comp.kj.params.old,
                                                      tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                      expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                      z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                      penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-burr"     = {temp = burr.params.m.recur(burr.params.old = comp.kj.params.old,
                                                      tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                      expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                      z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                      penalty = penalty, hyper.params = hyper.params.kj) },
          # Frequency distributions & their zero-inflation
          "poisson"     = {temp = poisson.params.m.recur(poisson.params.old = comp.kj.params.old,
                                                         tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                         expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                         z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                         penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-poisson"  = {temp = poisson.params.m.recur(poisson.params.old = comp.kj.params.old,
                                                         tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                         expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                         z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                         penalty = penalty, hyper.params = hyper.params.kj) },
          "nbinom"      = {temp = nbinom.params.m.recur(nbinom.params.old = comp.kj.params.old,
                                                        tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                        expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                        z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                        penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-nbinom"   = {temp = nbinom.params.m.recur(nbinom.params.old = comp.kj.params.old,
                                                        tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                        expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                        z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                        penalty = penalty, hyper.params = hyper.params.kj) },
          "binom"       = {temp = binom.params.m.recur(binom.params.old = comp.kj.params.old,
                                                        tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                        expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                        z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                        penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-binom"    = {temp = binom.params.m.recur(binom.params.old = comp.kj.params.old,
                                                        tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                        expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                        z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                        penalty = penalty, hyper.params = hyper.params.kj) },
          "gammacount"  = {temp = gammacount.params.m.recur(gammacount.params.old = comp.kj.params.old,
                                                            tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                            expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                            z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                            penalty = penalty, hyper.params = hyper.params.kj) },
          "ZI-gammacount"  = {temp = gammacount.params.m.recur(gammacount.params.old = comp.kj.params.old,
                                                               tl = tl.k, yl = yl.k, yu = yu.k, tu = tu.k,
                                                               expert.ll = comp.kj.pos.expert.ll, expert.tn = comp.kj.pos.expert.tn, expert.tn.bar = comp.kj.pos.expert.tn.bar,
                                                               z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e,
                                                               penalty = penalty, hyper.params = hyper.params.kj) },
          # Error
          stop("Invalid distribution!")
  )
  return(temp)
}
