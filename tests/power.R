library(elquantile)
near <- function(x, y, tol = 1e-7)
    stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
for (k in 1:4)
    near(elq_power(rep(100, k), rep(0, k), rep(.4, k))$power, .05)
n <- c(80, 120)
f <- c(.4, .2)
delta <- .5
variance <- .25 * sum(1 / (n * f^2))
expected <- pchisq(qchisq(.95, 1),
                   1,
                   ncp = delta^2 / variance,
                   lower.tail = FALSE)
near(elq_power(n, c(delta, 0), f)$power, expected)
near(elq_power(n, c(delta, 0), f, alternative = "greater")$power,
     pnorm(delta / sqrt(variance) - qnorm(.95)))
plan <- elq_sample_size(
    quantiles = c(0, .5),
    density = dnorm(0),
    ratio = c(1, 2)
)
stopifnot(plan$power >= .8, plan$previous.power < .8, plan$n[2] == 2 * plan$n[1])
for (m in c("el", "ael", "sign")) {
    z <- elq_power_exact(24, .3, method = m)
    near(z$power, sum(dbinom(0:24, 24, .3)[z$counts$reject]))
    near(z$null.size, sum(dbinom(0:24, 24, .5)[z$counts$reject]))
}
for (nn in c(11, 50))
    for (tt in c(.1, .37, .5, .9)) {
        pv <- elq_power_exact(nn, .3, tau = tt, method = "sign")$counts$p.value
        reference <- vapply(0:nn, function(mm)
            binom.test(mm, nn, tt)$p.value, numeric(1))
        near(pv, reference, 1e-12)
    }
plan.exact <- elq_sample_size_exact(
    power = .7,
    p = .3,
    method = "sign",
    max.n = 200
)
stopifnot(plan.exact$power >= .7, all(head(plan.exact$history$power, -1) < .7))
# A rare-quantile EL test can be liberal; exact enumeration must reveal it.
stopifnot(elq_power_exact(20, .1, tau = .1, method = "el")$null.size > .1)
set.seed(991)
saved <- .Random.seed
a <- elq_power_sim(30, function(n)
    rnorm(n, .4), nsim = 20, seed = 17)
b <- elq_power_sim(30, function(n)
    rnorm(n, .4), nsim = 20, seed = 17)
stopifnot(identical(saved, .Random.seed),
          identical(a$replicates, b$replicates))
broken <- elq_power_sim(20, function(n)
    rep(1, n), nsim = 4, seed = 1)
stopifnot(broken$failures == 4,
          is.na(broken$power),
          identical(unname(broken$failure.bounds), c(0, 1)))
paired <- function(n) {
    a <- rnorm(n)
    cbind(a, a + rnorm(n, .2))
}
sim <- elq_power_sim(25,
                     paired,
                     paired = TRUE,
                     nsim = 5,
                     seed = 3)
stopifnot(sim$valid == 5)
cat("Analytical, exact-count, simulation and sample-size checks passed.\n")
