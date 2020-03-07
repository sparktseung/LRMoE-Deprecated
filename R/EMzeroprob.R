## Computation of z
#' ECM: M-Step for \code{zero.prob}.
#'
#' @param z.zero.e.obs An object returned by \code{\link{z.zero.e.obs.recur}}.
#' @param z.pos.e.obs An object returned by \code{\link{comp.zkz.e.recur}}.
#' @param z.zero.e.lat An object returned by \code{\link{z.zero.e.lat.recur}}.
#' @param z.pos.e.lat An object returned by \code{\link{comp.zkz.e.recur}}.
#' @param k.e An object returned by \code{\link{comp.zkz.e.recur}}.
#'
#' @return \code{zero.prob} Updated zero.prob.
#'
#' @keywords internal
#'
# #' @export zero.prob.m.recur
zero.prob.m.recur = function(z.zero.e.obs, z.pos.e.obs, z.zero.e.lat, z.pos.e.lat, k.e)
{
  term.zero = sweep(matrix(z.zero.e.obs), 1, sweep(matrix(z.zero.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
    # z.zero.e.obs  + sweep(z.zero.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)
  term.pos  = sweep(matrix(z.pos.e.obs), 1, sweep(matrix(z.pos.e.lat), 1, matrix(k.e), FUN = "*", check.margin = FALSE), FUN = "+", check.margin = FALSE)
    # z.pos.e.obs   + sweep(z.pos.e.lat, 1, k.e, FUN = "*", check.margin = FALSE)

  numerator   = apply(term.zero, 2, sum)
  denominator = numerator + apply(term.pos, 2, sum) # apply(term.zero, 2, sum) + apply(term.pos, 2, sum)
  return( numerator / denominator )
}

