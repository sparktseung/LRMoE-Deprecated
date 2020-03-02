# Test NegativeBinomial
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
comp.dist = matrix(c("gammacount", "ZI-gammacount",
                     "ZI-gammacount", "gammacount"),
                   nrow = d, byrow = TRUE)
zero.prob = matrix(c(0, 0.3,
                     0.5, 0),
                   nrow = d, byrow = TRUE)
params.list = list(list(c(5, 0.2), c(10, 0.5)),
                   list(c(6, 0.3), c(8, 0.4)))

hyper.alpha = 5
hyper.params = list(list(c(2, 1, 2, 1), c(2, 1, 2, 1)),
                    list(c(2, 1, 2, 1), c(2, 1, 2, 1)))


Y.gc = dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)

head(Y.gc)

# Fit

alpha.guess = matrix(0, nrow = g, ncol = P)
comp.dist = comp.dist
zero.guess = matrix(c(0, 0.8,
                      0.3, 0),
                    nrow = d, byrow = TRUE)
params.guess = list(list(c(5, 0.5), c(9, 0.4)),
                    list(c(8, 0.2), c(7, 0.5)))

Y.obs = cbind(rep(0, N), Y.gc[,1], Y.gc[,1], rep(Inf, N),
              rep(0, N), Y.gc[,2], Y.gc[,2], rep(Inf, N))


model.fit = LRMoE.fit(Y = Y.obs, X = X,
                      n.comp = 2, comp.dist = comp.dist,
                      alpha.init = alpha.guess,
                      zero.init = zero.guess,
                      params.init = params.guess,
                      penalty = TRUE, hyper.alpha = hyper.alpha, hyper.params = hyper.params)






