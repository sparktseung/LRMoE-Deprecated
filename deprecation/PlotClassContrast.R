# #' Plots two stacked bar charts of latent class proportions, which contrasts prior and posterior latent class proportions, given a covariate matrix.
# #'
# #' @param X A matrix of covariates.
# #' @param alpha.old,alpha.new A matrix of logit regression coeffiecients.
# #'
# #' @return A \code{ggplot2} object.
# #'
# #' @seealso \code{\link{LRMoE.fit}}, \code{\link{predict.class}}
# #'
# #' @import ggplot2
# #'
# #' @export plot.dataset.contrast.prob
# #'

# plot.dataset.contrast.prob = function(X, alpha.old, alpha.new)
# {
#   # Old assignments
#   assign.class.old = predict.class(X, alpha.old)
#   df.old = data.frame(class = assign.class.old, count = rep(1, length(assign.class.old)))
#   df.old$class = factor(df.old$class, levels = c(1:nrow(alpha.old)))
#
#   df.old.aggre = aggregate(df.old$count, by = list(class = df.old$class), FUN = "sum")
#   df.old.plot = data.frame(cbind(case = rep("Prior", nrow(df.old.aggre)), df.old.aggre) )
#   df.old.plot$probability = df.old.plot[,3]/sum(df.old.plot[,3])
#
#   # New assignments
#   assign.class.new = predict.class(X, alpha.new)
#   df.new = data.frame(class = assign.class.new, count = rep(1, length(assign.class.new)))
#   df.new$class = factor(df.new$class, levels = c(1:nrow(alpha.old)))
#
#   df.new.aggre = aggregate(df.new$count, by = list(class = df.new$class), FUN = "sum")
#   df.new.plot = data.frame(cbind(case = rep("Posterior", nrow(df.new.aggre)), df.new.aggre) )
#   df.new.plot$probability = df.new.plot[,3]/sum(df.new.plot[,3])
#
#   # df to plot
#   df.plot = data.frame(rbind(df.old.plot, df.new.plot))
#   df.plot$case = factor(df.plot$case, levels = c("Prior", "Posterior"))
#
#   return(
#     ggplot(df.plot, aes(fill=class, y=probability, x = case)) +
#       geom_bar(position="fill", stat="identity") +
#       # geom_text(aes(label = round(probability, 2)),
#       #           position = position_stack(vjust = 0.5)) +
#       xlab("") + ylab("Latent Class Proportion") +
#       # ggtitle("Proportion of Latent Classes") +
#       theme(plot.title = element_text(hjust = 0.5))
#   )
# }
