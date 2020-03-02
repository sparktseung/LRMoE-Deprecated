# Test the dataset simulator
N = 100
P = 5
g = 3
d = 2

X = matrix(rnorm(N*P), nrow = N)
alpha = matrix(rnorm(g*P), nrow = g)

comp.dist = matrix(c("gamma", "ZI-lnorm", "invgauss",
                     "ZI-poisson", "nbinom", "ZI-gammacount"),
                   nrow = d, byrow = TRUE)

zero.prob = matrix(c(0, 0.4, 0,
                     0.3, 0, 0.8),
                   nrow = d, byrow = TRUE)

params.list = list(list(c(3, 0.5), c(0.3, 0.2), c(5, 2)),
                   list(c(5), c(2, 0.5), c(3, 0.5))
                   )

Y = dataset.simulator(X, alpha, comp.dist, zero.prob, params.list)

summary(Y)
hist(Y[,1])
hist(Y[,2])


# Test loglikelihood
Y.obs = cbind(rep(0,N), Y[,1], Y[,1], rep(Inf, N), rep(0,N), Y[,2], Y[,2], rep(Inf, N))
hyper.alpha = 5
hyper.params = list(list(c(2, 5, 2, 5), c(2, 5, 2, 5), c(2, 5, 2, 5)),
                    list(c(2, 5), c(2, 5), c(2, 5, 2, 5)))
gate.ll = gate.logit(X, alpha)
expert.list = expert.loglik.dim.comp(Y.obs, comp.dist, zero.prob, params.list)
ll.list = gate.expert.loglik(alpha, gate.ll, expert.list, penalty = TRUE, hyper.alpha, hyper.params)

model = list(alpha = alpha, comp.dist = comp.dist, zero.prob = zero.prob, params.list = params.list)

temp = LRMoE.loglik(X, Y.obs, model, penalty = TRUE, hyper.alpha, hyper.params)


