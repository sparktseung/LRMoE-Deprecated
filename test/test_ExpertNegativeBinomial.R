# Test the expert function outputs the correct values

N = 100

tl = rep(c(0, 1), 4*N)
yl = rep(c(1, 2, 3, 4), 2*N)
yu = rep(c(1, 2, 4, 4), 2*N)
tu = rep(c(4, 5, Inf, Inf), 2*N)
g = 1
size.n = 5
prob.p = 0.2

temp = expert.nbinom(tl, yl, yu, tu, g = 1, size.n, prob.p)

head(cbind(temp[[1]], temp[[2]], temp[[3]]), 8)
