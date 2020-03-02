# Test Weibull
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
comp.dist = matrix(c("weibull", "ZI-weibull",
                     "ZI-weibull", "weibull"),
                   nrow = d, byrow = TRUE)
zero.prob = matrix(c(0, 0.3,
                     0.5, 0),
                   nrow = d, byrow = TRUE)
params.list = list(list(c(2, 5), c(1, 10)),
                   list(c(3, 6), c(2, 15)))

hyper.alpha = 5
hyper.params = list(list(c(2, 1, 2, 1), c(2, 1, 2, 1)),
                    list(c(2, 1, 2, 1), c(2, 1, 2, 1)))


Y.weibull = dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)

head(Y.weibull)

# Fit

alpha.guess = matrix(0, nrow = g, ncol = P)
comp.dist = comp.dist
zero.guess = matrix(c(0, 0.8,
                      0.3, 0),
                    nrow = d, byrow = TRUE)
params.guess = list(list(c(1.5, 6), c(2.5, 5)),
                    list(c(2, 3), c(5, 20)))

Y.obs = cbind(rep(0, N), Y.weibull[,1], Y.weibull[,1], rep(Inf, N),
              rep(0, N), Y.weibull[,2], Y.weibull[,2], rep(Inf, N))


model.fit = LRMoE.fit(Y = Y.obs, X = X,
                      n.comp = 2, comp.dist = comp.dist,
                      alpha.init = alpha.guess,
                      zero.init = zero.guess,
                      params.init = params.guess,
                      penalty = TRUE, hyper.alpha = hyper.alpha, hyper.params = hyper.params)






