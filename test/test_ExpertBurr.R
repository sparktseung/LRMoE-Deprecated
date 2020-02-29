# Test the expert function outputs the correct values

N = 100

tl = rep(c(0, 1), 4*N)
yl = rep(c(1, 2, 3, 4), 2*N)
yu = rep(c(1, 2, 4, 4), 2*N)
tu = rep(c(4, 5, Inf, Inf), 2*N)
g = 1
shape1.k = 2
shape2.c = 4
scale.lambda = 2

temp = expert.burr(tl, yl, yu, tu, g = 1, shape1.k, shape2.c, scale.lambda)

head(cbind(temp[[1]], temp[[2]], temp[[3]]), 8)
