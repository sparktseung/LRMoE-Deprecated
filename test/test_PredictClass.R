# Test Predict.Class
library(LRMoE)
data("DemoData")
set.seed(777)

# Generate data
P = 5
g = 2
d = 2
N = 10000

X[, "agedriver"] = (X[, "agedriver"] - mean(X[, "agedriver"]))/sd(X[, "agedriver"])
X[, "agecar"] = (X[, "agecar"] - mean(X[, "agecar"]))/sd(X[, "agecar"])

alpha = matrix(rnorm(g*P), nrow = g)
comp.dist = matrix(c("lnorm", "ZI-lnorm",
                     "ZI-lnorm", "lnorm"),
                   nrow = d, byrow = TRUE)
zero.prob = matrix(c(0, 0.3,
                     0.5, 0),
                   nrow = d, byrow = TRUE)
params.list = list(list(c(2, 0.6), c(1, 0.5)),
                   list(c(1, 0.5), c(2, 0.4)))

hyper.alpha = 5
hyper.params = list(list(c(2, 1, 2, 1), c(2, 1, 2, 1)),
                    list(c(2, 1, 2, 1), c(2, 1, 2, 1)))


Y.lnorm = dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)

head(Y.gc)

# Fit

alpha.guess = matrix(0, nrow = g, ncol = P)
comp.dist = comp.dist
zero.guess = matrix(c(0, 0.8,
                      0.3, 0),
                    nrow = d, byrow = TRUE)
params.guess = list(list(c(1, 0.5), c(2, 0.4)),
                    list(c(2, 0.2), c(1, 0.5)))

Y.obs = cbind(rep(0, N), Y.lnorm[,1], Y.lnorm[,1], rep(Inf, N),
              rep(0, N), Y.lnorm[,2], Y.lnorm[,2], rep(Inf, N))


model.fit = LRMoE.fit(Y = Y.obs, X = X,
                      n.comp = 2, comp.dist = comp.dist,
                      alpha.init = alpha.guess,
                      zero.init = zero.guess,
                      params.init = params.guess,
                      penalty = TRUE, hyper.alpha = hyper.alpha, hyper.params = hyper.params)

# predict.class.prob
predict.class.prob(head(X), model.fit$alpha.fit)

plot.ind.class.prob(X[1,], model.fit$alpha.fit) # + ggtitle("haha")

# predict.class
predict.class(head(X), model.fit$alpha.fit)

plot.dataset.class.prob(X, model.fit$alpha.fit)

# predict.class.prob.posterior
predict.class.prob.posterior(head(X), head(Y.obs),
      model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)

plot.ind.class.prob.posterior(X[1,], Y.obs[1,],
      model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)

# predict.class.posterior
predict.class.posterior(head(X), head(Y.obs),
      model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)

plot.dataset.prob.posterior(X, Y.obs,
      model.fit$alpha.fit, model.fit$comp.dist, model.fit$zero.fit, model.fit$params.fit)




