#' Plot the fitted density, given a vector of covariates and a model
#'
#' @param X A vector of covariates.
#' @param alpha A g*P matrix of numerics, which contains initial guess of the logit regression coefficients.
#'                   The last row should all be zero, representing the default latent class.
#' @param comp.dist A d*g matrix of strings, which specify the component distributions to fit.
#'                  The rows represent the dimensions of \code{Y}, while the columns represent the component distributions.
#' @param zero.prob A d*g matrix of numerics in (0,1), which specify the probability mass at zero for component distributions.
#' @param params.list A list of length d, where each element is a sublist of length g.
#'                    Each sublist contains one numeric vector, which is the parameter value for the corresponding \code{comp.dist}.
#' @param plot.dim A numeric indicating which dimention of y to plot.
#' @param plot.lim Upper bound of y for plotting. Default is 50, if no value is provided.
#'
#' @return A \code{\link{ggplot2}} object.
#'
#' @seealso \code{\link{LRMoE.fit}}
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @rawNamespace S3method(plot, ind.fitted.dist)
#'
#' @export plot.ind.fitted.dist
#'
plot.ind.fitted.dist = function(X, alpha, comp.dist, zero.prob, params.list, plot.dim = 1, plot.lim = NULL)
{
  freq.dist = c("poisson", "ZI-poisson", "nbinom", "ZI-nbinom", "gammacount", "ZI-gammacount")

  weighting = predict.class.prob(X, alpha)

  # dim.m = nrow(comp.dist)
  dim.m = 1
  n.comp = ncol(comp.dist)

  if(is.null(plot.lim)){
    # plot.lim = rep(50, dim.m)
    plot.lim = 50
  }

  # result = list()

  for(k in plot.dim:plot.dim)
    # for(k in 2:3)
  {
    # Zero inflation
    zero.prob.d = as.numeric( zero.prob[k,] %*% t(weighting) )

    # Positive density series
    y.series = NULL
    if(comp.dist[k,1] %in% freq.dist){ # Count distributions: 1
      # y.series = seq(from = 0, to = plot.lim[k], by = 1)
      y.series = seq(from = 0, to = plot.lim, by = 1)
    }else{ # Frequency distributions: 0.05
      # y.series = seq(from = 0, to = plot.lim[k], by = 0.05)
      y.series = seq(from = 0, to = plot.lim, by = 0.05)
    }

    dens.series = matrix(0, ncol = n.comp, nrow = length(y.series))
    for(j in 1:n.comp){
      dens.series[,j] = ind.dens.y.pos(comp.dist[k,j], params.list[[k]][[j]], y.series)
    }

    # Some adjustment of zero.prob.d for count distributions
    pos.dens.series = dens.series %*% t(weighting)
    # all.dens.series = pos.dens.series * (1-zero.prob.d)
    all.dens.series = pos.dens.series
    # all.dens.series[1] = all.dens.series[1] + zero.prob.d

    # Plotting
    if(comp.dist[k,1] %in% freq.dist){ # Count distributions: move the zero mass
      # zero.prob.d = zero.prob.d + pos.dens.series[1,1]

      # Getting together the data.frame
      df = data.frame(y = y.series)
      # df$label = rep("pos", nrow(df))
      # df$label[1] = "zero"
      # df$label = factor(df$label, levels = c("zero", "pos"))
      df$zero = c(zero.prob.d, rep(0, nrow(df)-1))
      df$pos = c(all.dens.series)

      df1 = melt(df, id.var = "y")
      colnames(df1) = c("y.series", "label", "prob")

      temp =
        ggplot(df1, aes(x = y.series, y = prob, fill = label)) +
        geom_bar(stat='identity', position = "stack", width = 0.75) +
        xlab("Y") + ylab("Probability Mass") +
        ggtitle(paste0("Fitted Distribution of Dimension ", k)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(legend.position='none')

    }else{
      # Getting together the data.frame

      df = data.frame(y = y.series)
      df$label = rep("pos", nrow(df))
      df$label[1] = "zero"
      df$label = factor(df$label, levels = c("zero", "pos"))
      df$zero = c(zero.prob.d, rep(0, nrow(df)-1))
      df$pos = c(0, all.dens.series[-1])

      temp =
        ggplot(df, aes(x = y.series)) +
        geom_bar(stat = "identity", aes(y = zero), fill = c("#00BFC4"), width = 0.01*plot.lim, position = "identity") +
        xlab("Y") + ylab("Zero Inflation") +
        geom_line(aes(y = pos/max(df$pos)), color = "#F8766D", size = 1) +
        scale_y_continuous(sec.axis = sec_axis(~.*max(df$pos), name = "Positive Density")) +
        ggtitle(paste0("Fitted Distribution of Dimension ", k)) +
        theme(plot.title = element_text(hjust = 0.5))
    }

  }

  return(temp)
}
