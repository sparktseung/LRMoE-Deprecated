## Computation of z, k, z
#' ECM: E-Step for \code{z.e.obs}, \code{k.e} and \code{z.e.lat}.
#'
#' @param gate.ll An object returned by \code{\link{gate.logit}}.
#' @param expert.list An object returned by \code{\link{expert.loglik.dim.comp}}.
#' @param ll.list An object returned by \code{\link{gate.expert.loglik}}.
#'
#' @return \code{z.e.obs},\code{k.e},\code{z.e.lat} Numerical vectors of length N.
#'
#' @keywords internal
#'
#' @importFrom matrixStats rowLogSumExps
#'
#' @export comp.zkz.e.recur
comp.zkz.e.recur = function(gate.ll, expert.list, ll.list)
{

  sample.size.n = ll.list$sample.size.n # sample size
  n.comp = ll.list$n.comp # no. of experts
  dim.m = ll.list$dim.m # no. of dimensions of y

  # Note: For the purpose of optimizing wrt gating weights alpha's,
  #       there is no need to consider z.zero and z.pos yet.
  z.e.obs       = array(0, dim = c(sample.size.n, n.comp))
  # z.zero.obs  = array(0, dim = c(sample.size.n, n.comp))
  # z.pos.obs   = array(0, dim = c(sample.size.n, n.comp))

  z.e.lat       =  rep(list(array(0, dim = c(sample.size.n, n.comp))), dim.m)
  # z.zero.lat  = array(0, dim = c(sample.size.n, n.comp))
  # z.pos.lat   = array(0, dim = c(sample.size.n, n.comp))

  k.e           = array(0, dim = c(sample.size.n, 1))

  z.e.obs = exp( sweep(ll.list$ll.ind, 1, ll.list$comp.aggre.ll.ind, FUN = "-", check.margin = FALSE) )

  z.e.lat = exp( sweep(ll.list$ll.ind.tn.bar, 1, ll.list$comp.aggre.ll.ind.tn.bar, FUN = "-", check.margin = FALSE) )
  z.e.lat[is.na(z.e.lat)] = 1/n.comp

  return(list(z.e.obs = z.e.obs, z.e.lat = z.e.lat, k.e = k.e))

}
