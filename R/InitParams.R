## Severity Initialization
#' Initializes parameter for severity distributions using CMM.
#'
#' @param Y A vector of response variables.
#' @param cluster The \code{cluster} list vector returned by \code{\link[stats]{kmeans}}
#' @return A list of parameter initialization.
#'
#' @importFrom EnvStats skewness kurtosis
#' @importFrom stats var
#'
#' @export cluster.mm.severity
cluster.mm.severity = function(Y, cluster)
{
  n.group = length(unique(cluster))

  result = NULL

  for(j in 1:n.group){
    subset = Y[which(cluster==j)]
    subset.pos = subset[which(subset>0)]

    cluster.prop = length(subset) / length(Y)

    zero.prop = 1 - length(subset.pos)/length(subset)

    mean.pos = mean(subset.pos)
    var.pos = var(subset.pos)
    cv.pos = sqrt(var.pos)/mean.pos
    skew.pos = skewness(subset.pos)
    kurt.pos = kurtosis(subset.pos)

    gamma.init = c(shape = mean.pos^2/var.pos, scale = var.pos/mean.pos)
    lnorm.init = c(meanlog = log(mean.pos)- 0.5*(log(var.pos/(mean.pos^2) + 1)), sdlog = sqrt(log(var.pos/(mean.pos^2) + 1)))
    invgauss.init = c(mean = mean.pos, scale = (mean.pos)^3/var.pos )
    weibull.init = c(shape = 3, scale = mean.pos * 1.5 )# Ad-hoc
    burr.init = c(shape1 = 2, shape2 = 5, scale = mean.pos/0.136) # Ad-hoc. 0.136 = Beta(4.8, 1.2)

    result[[length(result)+1]] = list(cluster.prop = cluster.prop, zero.prop = zero.prop, mean.pos = mean.pos, var.pos = var.pos, cv.pos = cv.pos, skew.pos = skew.pos, kurt.pos = kurt.pos,
                                      gamma.init = gamma.init, lnorm.init = lnorm.init, invgauss.init = invgauss.init,
                                      weibull.init = weibull.init, burr.init = burr.init)
  }

  return(result)
}


## Frequency Initialization
#' Initializes parameter for frequency distributions using CMM.
#'
#' @param Y A vector of response variables.
#' @param cluster The \code{cluster} list vector returned by \code{\link[stats]{kmeans}}
#' @return A list of parameter initialization.
#'
#' @importFrom EnvStats skewness kurtosis
#' @importFrom stats var
#'
#' @export cluster.mm.frequency
cluster.mm.frequency = function(Y, cluster)
{
  n.group = length(unique(cluster))

  result = NULL

  for(j in 1:n.group){
    subset = Y[which(cluster==j)]
    subset.pos = subset[which(subset>0)]

    cluster.prop = length(subset) / length(Y)

    zero.prop = 1 - length(subset.pos)/length(subset)

    mean.pos = mean(subset.pos)
    var.pos = var(subset.pos)
    cv.pos = sqrt(var.pos)/mean.pos
    skew.pos = skewness(subset.pos)
    kurt.pos = kurtosis(subset.pos)

    poisson.init = c(lambda = mean.pos)
    nbinom.init = c(size.n = mean.pos * (mean.pos/var.pos/(1-mean.pos/var.pos)), prob.p = mean.pos/var.pos)
    gammacount.init = c(shape = mean.pos * (mean.pos/var.pos/(1-mean.pos/var.pos)), scale = mean.pos/var.pos ) # ad-hoc, same as nbinom


    result[[length(result)+1]] = list(cluster.prop = cluster.prop, zero.prop = zero.prop, mean.pos = mean.pos, var.pos = var.pos, cv.pos = cv.pos, skew.pos = skew.pos, kurt.pos = kurt.pos,
                                      poisson.init = poisson.init, nbinom.init = nbinom.init, gammacount.init = gammacount.init)
  }

  return(result)
}


