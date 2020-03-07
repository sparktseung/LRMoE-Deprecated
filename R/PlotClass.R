#' Plots a stacked bar chart of most likely latent class proportion, given a matrix of covariates.
#'
#' @param X A matrix of covariates.
#' @param alpha A matrix of logit regression coeffiecients.
#' @param title A text string for plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{LRMoE.fit}}, \code{\link{predict.class}}
#'
#' @import ggplot2
#' @importFrom stats aggregate
#'
#' @rawNamespace S3method(plot, dataset.class.prob)
#'
#' @export plot.dataset.class.prob
#'
plot.dataset.class.prob = function(X, alpha, title = "Proportion of Latent Classes")
{
  assign.class = predict.class(X, alpha)
  df = data.frame(class = assign.class, count = rep(1, length(assign.class)))
  df$class = factor(df$class, levels = c(1:nrow(alpha)))

  df.aggre = aggregate(df$count, by = list(class = df$class), FUN = "sum")
  df.plot = data.frame(cbind(covariate = rep("X", nrow(df.aggre)), df.aggre) )
  df.plot$probability = df.plot[,3]/sum(df.plot[,3])

  return(
    ggplot(df.plot, aes(fill=class, y=probability, x = covariate)) +
      geom_bar(position="fill", stat="identity") +
      # geom_text(aes(label = round(probability, 2)),
      #           position = position_stack(vjust = 0.5)) +
      xlab("") + ylab("Latent Class Proportions") +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}
