# Test the gating function outputs the correct values
library(matrixStats)

N = 10 # Sample size
P = 5 # Number of covariates
g = 3 # Number of latent classes

x = matrix(rnorm(N*P), nrow = N, ncol = P)
alpha = matrix(rnorm(g*P), nrow = g, ncol = P)

temp = gate.logit(x, alpha)

# Check dimension
dim(temp)

# Check row sum
exp(rowLogSumExps(temp))


