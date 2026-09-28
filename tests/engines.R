library(elquantile)
near <- function(x, y, tol = 1e-7)
    stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
for (tau in c(.1, .5, .9))
    for (m in c(1, 4, 10, 17, 19)) {
        n <- 20
        fit <- elq_dual(c(1 - tau, -tau), c(m, n - m))
        expected <- 2 * (m * log((m / n) / tau) + (n - m) * log((1 - m / n) /
                                                                    (1 - tau)))
        near(fit$statistic, expected, 1e-7)
        near(sum(fit$probabilities), 1)
        stopifnot(abs(fit$residual) < 1e-7)
    }
stopifnot(!elq_dual(rep(1, 5))$feasible, is.infinite(elq_dual(rep(1, 5))$statistic))
for (n in c(2, 10, 100)) {
    a <- max(1, log(n) / 2)
    expected <- -2 * (n * log((n + 1) * a / (n * (a + 1))) + log((n + 1) /
                                                                     (a + 1)))
    fit <- elq_dual(rep(.4, n), adjustment = a)
    near(fit$statistic, expected)
    near(sum(fit$probabilities), 1)
    stopifnot(abs(fit$residual) < 1e-8)
}
for (kernel in c("biweight", "epanechnikov", "zhou_jing")) {
    near(integrate(function(u)
        elq_kernel(u, kernel), -1, 1)$value, 1)
    near(integrate(function(u)
        u * elq_kernel(u, kernel), -1, 1)$value, 0)
    near(elq_kernel(c(-2, 2), kernel, TRUE), c(0, 1))
}
# This special cancellation distinguishes Zhou-Jing's kernel from standard KDE.
near(integrate(
    function(u)
        u * elq_kernel(u, "zhou_jing") *
        elq_kernel(u, "zhou_jing", TRUE),
    -1,
    1
)$value, 0)
stopifnot(elq_kernel(.99, "zhou_jing") < 0, max(elq_kernel(
    seq(-1, 1, length.out = 1001), "zhou_jing", TRUE
)) > 1)

x <- c(-2, -.5, .1, 1, 3)
theta <- c(-2, -1, .1, 2.9)
p <- approx(x, (1:5 - .5) / 5, theta)$y
expected <- 10 * (p * log(p / .5) + (1 - p) * log((1 - p) / .5))
near(elq_statistic(x, theta, method = "adimari"), expected)
stopifnot(is.infinite(elq_statistic(x, 3, method = "adimari")))
set.seed(813)
x <- rnorm(61)
h <- .3
q <- .2
b <- elq_bartlett(x, bandwidth = h)
raw <- elq_statistic(x, q, method = "chen_hall", bandwidth = h)
corrected <- elq_statistic(x, q, method = "chen_hall_bartlett", bandwidth =
                               h)
near(corrected, raw / (1 + b / length(x)))
near(elq_bartlett(x, type = "limit"), .5)
# Affine equivariance checks bandwidths, interpolation, and parameter units.
for (m in elq_methods()$method) {
    w <- elq_statistic(x, q, method = m)
    w2 <- elq_statistic(3 * x + 7, 3 * q + 7, method = m)
    near(w, w2, 1e-6)
}
cat("Engine identities and kernel checks passed.\n")
