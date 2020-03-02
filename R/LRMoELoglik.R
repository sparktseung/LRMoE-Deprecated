## Computes the loglikelihood
#' Computes the loglikelihood of LRMoE, given X, Y and a fitted model.
#'
#' @param X A N*P matrix of covariates.
#' @param Y A N*d matrix of response.
#' @param model A list of parameters specifying an LRMoE model, including
#' \itemize{
#'      \item \code{alpha}: A g*P matrix, where g is the number of components and P is the number of covariates.
#'      \item \code{comp.dist}: A d*g matrix of strings, describing component distributions by dimension and by component.
#'      \item \code{zero.prob}: A d*g matrix of numbers between 0 and 1, describing zero probability masses by dimension and by component.
#'      \item \code{params.list}: A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the initial parameter guess for the corresponding \code{comp.dist}.
#'        }
#' @param penalty TRUE/FALSE, which indicates whether parameter penalty should be applied. Default (and recommended) is TRUE.
#' @param hyper.alpha A numeric, which penalizes the magnitude of \code{alpha}.
#' @param hyper.params A list of length d. Each element is a sublist of length g.
#'                     Each element of a sublist is a vector of numerics, which penalizes expert parameters. See also \code{\link{expert.loglik.pen.dim.comp}}.
#'
#' @return Loglikelihood (with and without penalty), AIC and BIC.
#'
#' @export LRMoE.loglik
#'
LRMoE.loglik = function(X, Y, model, penalty = TRUE, hyper.alpha, hyper.params)
{
  sample.size.n = nrow(X)
  gate.ll = gate.logit(X, model$alpha)
  expert.list = expert.loglik.dim.comp(Y, model$comp.dist, model$zero.prob, model$params.list)
  ll.list = gate.expert.loglik(model$alpha, gate.ll, expert.list, penalty = TRUE, hyper.alpha, hyper.params)

  AIC = -2*ll.list$ll.np + 2*(count.alpha(model$alpha) + count.zero(model$comp.dist) + count.pos.params(model$comp.dist))
  BIC = -2*ll.list$ll.np + log(sample.size.n)*(count.alpha(model$alpha) + count.zero(model$comp.dist) + count.pos.params(model$comp.dist))

  return(list(ll = ll.list$ll, ll.np = ll.list$ll.np, AIC = AIC, BIC = BIC))
}
