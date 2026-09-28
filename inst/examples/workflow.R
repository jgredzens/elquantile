library(elquantile)
set.seed(27)
x <- rexp(60)
y <- rexp(80) + .3
print(elq_one_sample(
    x,
    q0 = log(2),
    method = "chen_hall_bartlett",
    conf.int = TRUE
))
print(elq_two_sample(x, y, method = "ael"))
print(elq_posthoc(list(
    A = x, B = y, C = rexp(70)
), method = "el"))
print(elq_sample_size(
    quantiles = c(0, .5),
    density = c(dnorm(0), dnorm(0) / 2),
    ratio = c(1, 2)
))
print(elq_power_exact(40, pnorm(0, mean = .5), method = "el"))
print(elq_power_sim(40, function(n)
    rnorm(n, .5), method = "adimari", nsim = 100, seed = 29))
