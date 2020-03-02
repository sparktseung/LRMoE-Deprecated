# #' Plots two stacked bar charts of latent class probabilities, which contrasts prior and posterior latent class probabilities, given a covariate vector.
# #'
# #' @param X A vector of covariates.
# #' @param alpha.old,alpha.new A matrix of logit regression coeffiecients.
# #'
# #' @return A \code{ggplot2} object.
# #'
# #' @seealso \code{\link{LRMoE.fit}}, \code{\link{predict.class.prob}}
# #'
# #' @import ggplot2
# #'
# #' @export plot.ind.prob.contrast
# #'

# plot.ind.prob.contrast = function(X, alpha.old, alpha.new)
# {
#   weighting.old = predict.class.prob(X, alpha.old)
#   weighting.new = predict.class.prob(X, alpha.new)
#
#   df = data.frame(case = c(rep("Prior", nrow(alpha.old)), rep("Posterior", nrow(alpha.old))), class = rep(as.factor(c(1:nrow(alpha.old))), 2), probability = (c(weighting.old, weighting.new)))
#   df$case = factor(df$case, levels = c("Prior", "Posterior"))
#
#   return(
#     ggplot(df, aes(fill=class, y=probability, x = case)) +
#       geom_bar(position="fill", stat="identity") +
#       # geom_text(aes(label = round(probability, 2)),
#       #           position = position_stack(vjust = 0.5)) +
#       xlab("") + ylab("Latent Class Probabilities") +
#       # ggtitle("Prediction of Latent Classes") +
#       theme(plot.title = element_text(hjust = 0.5))
#   )
# }
