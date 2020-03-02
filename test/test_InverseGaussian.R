# Test InverseGaussian
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
comp.dist = matrix(c("invgauss", "ZI-invgauss",
                     "ZI-invgauss", "invgauss"),
                   nrow = d, byrow = TRUE)
zero.prob = matrix(c(0, 0.3,
                     0.5, 0),
                   nrow = d, byrow = TRUE)
params.list = list(list(c(10, 2), c(5, 1)),
                   list(c(6, 5), c(9, 4)))

hyper.alpha = 5
hyper.params = list(list(c(2, 1, 2, 1), c(2, 1, 2, 1)),
                    list(c(2, 1, 2, 1), c(2, 1, 2, 1)))


Y.invgauss = dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)

head(Y.invgauss)

# Fit

alpha.guess = matrix(0, nrow = g, ncol = P)
comp.dist = comp.dist
zero.guess = matrix(c(0, 0.8,
                      0.3, 0),
                    nrow = d, byrow = TRUE)
params.guess = list(list(c(8, 1), c(6, 2)),
                    list(c(5, 2), c(10, 5)))

Y.obs = cbind(rep(0, N), Y.invgauss[,1], Y.invgauss[,1], rep(Inf, N),
              rep(0, N), Y.invgauss[,2], Y.invgauss[,2], rep(Inf, N))


model.fit = LRMoE.fit(Y = Y.obs, X = X,
                      n.comp = 2, comp.dist = comp.dist,
                      alpha.init = alpha.guess,
                      zero.init = zero.guess,
                      params.init = params.guess,
                      penalty = TRUE, hyper.alpha = hyper.alpha, hyper.params = hyper.params)






