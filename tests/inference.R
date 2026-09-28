library(elquantile)
near <- function(x, y, tol = 1e-6)
    stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
set.seed(317)
x <- rnorm(35)
y <- rnorm(35) + .3 * x
z <- rnorm(40, sd = 2)
ctrl <- elq_control(grid.points = 41, warn = FALSE)
for (m in elq_methods()$method) {
    one <- elq_one_sample(x,
                          q0 = .1,
                          method = m,
                          control = ctrl)
    paired <- elq_two_sample(
        x,
        y,
        delta = .1,
        paired = TRUE,
        method = m,
        control = ctrl
    )
    direct <- elq_one_sample(x - y,
                             q0 = .1,
                             method = m,
                             control = ctrl)
    near(paired$statistic, direct$statistic)
    two <- elq_two_sample(x, y, method = m, control = ctrl)
    grouped <- elq_k_sample(list(x = x, y = y), method = m, control = ctrl)
    near(two$statistic, grouped$statistic)
    k <- elq_k_sample(list(x = x, y = y, z = z), method = m, control = ctrl)
    stopifnot(
        one$p.value >= 0,
        one$p.value <= 1,
        k$parameter == 2,
        k$statistic >= 0,
        two$statistic >= 0
    )
    ci <- elq_confint(x, method = m, control = ctrl)
    candidates <- seq(min(x), max(x), length.out = 201)
    accepted <- elq_statistic(x, candidates, method = m, control = ctrl) <= qchisq(.95, 1)
    stopifnot(all(candidates[accepted] >= ci[1] - 1e-6), all(candidates[accepted] <= ci[2] +
                                                                 1e-6))
    if (m %in% c("chen_hall", "chen_hall_bartlett")) {
        near(elq_statistic(
            x,
            as.numeric(ci),
            method = m,
            control = ctrl
        ),
        rep(qchisq(.95, 1), 2),
        1e-5)
    }
}
for (m in c("el", "ael")) {
    ci <- elq_confint(x,
                      y,
                      method = m,
                      interval = "profile",
                      control = ctrl)
    deltas <- seq(ci[1] + .01, ci[2] - .01, length.out = 10)
    stopifnot(all(
        vapply(deltas, function(d)
            elq_two_sample(
                x,
                y,
                delta = d,
                method = m,
                control = ctrl
            )$p.value >= .05, logical(1))
    ))
}
ph <- elq_posthoc(
    list(A = x, B = y, C = z),
    tau = c(.25, .5, .75),
    conf.int = TRUE,
    control = ctrl
)
stopifnot(nrow(ph) == 9, attr(ph, "family.size") == 9)
near(ph$p.adj, p.adjust(ph$p.value, "holm"))
near(elq_sign_test(x, q0 = .1)$p.value,
     binom.test(sum(x <= .1), length(x), .5)$p.value)
stopifnot(elq_sign_test(x + 5, alternative = "greater")$p.value < .001)
fails <- function(expr)
    inherits(tryCatch(
        force(expr),
        error = function(e)
            e), "error")
stopifnot(fails(elq_test(c(1, 1, 2))),
          fails(elq_test(c(1, NA, 2))),
          fails(elq_test(list(x, y, z), alternative = "greater")),
          fails(elq_two_sample(x, y[-1], paired = TRUE)))
cat("Inference, pairing, confidence sets and multiplicity checks passed.\n")
