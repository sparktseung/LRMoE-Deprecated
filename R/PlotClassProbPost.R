#' Plots two stacked bar charts of latent class probabilities, which contrasts prior and posterior latent class probabilities, given a covariate vector.
#'
#' @param X A matrix of covariates.
#' @param Y A matrix of observed responses for \code{X}.
#' @param alpha A matrix of logit regression coefficients.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{LRMoE.fit}}, \code{\link{predict.class.prob.posterior}}
#'
#' @import ggplot2
#'
#' @export plot.ind.class.prob.posterior
#'
plot.ind.class.prob.posterior = function(X, Y, alpha, comp.dist, zero.prob, params.list, title = "Prediction of Latent Classes")
{
  # prior
  alpha.old = alpha
  weighting.old = predict.class.prob(X, alpha.old)

  # posterior
  Y = matrix(Y, nrow = 1)
  weighting.new = predict.class.prob.posterior(X, Y, alpha, comp.dist, zero.prob, params.list)

  df = data.frame(case = c(rep("Prior", nrow(alpha.old)), rep("Posterior", nrow(alpha.old))), class = rep(as.factor(c(1:nrow(alpha.old))), 2), probability = (c(weighting.old, weighting.new)))
  df$case = factor(df$case, levels = c("Prior", "Posterior"))

  return(
    ggplot(df, aes(fill=class, y=probability, x = case)) +
      geom_bar(position="fill", stat="identity") +
      # geom_text(aes(label = round(probability, 2)),
      #           position = position_stack(vjust = 0.5)) +
      xlab("") + ylab("Latent Class Probabilities") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))

  )
}
