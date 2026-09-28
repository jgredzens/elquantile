# All values are -2 log(relative likelihood).
.eq_kl <- function(p, tau) {
    ans <- rep(Inf, length(p))
    ok <- is.finite(p) & p >= 0 & p <= 1
    z <- p[ok]
    a <- b <- numeric(length(z))
    positive <- z > 0
    below <- z < 1
    a[positive] <- z[positive] * log(z[positive] / tau)
    b[below] <- (1 - z[below]) * log((1 - z[below]) / (1 - tau))
    ans[ok] <- pmax(0, a + b)
    ans
}
.eq_binary <- function(m, n, tau) {
    value <- 2 * n * .eq_kl(m / n, tau)
    value[m == 0 | m == n] <- Inf
    value
}

elq_dual <- function(g,
                     weights = rep(1, length(g)),
                     adjustment = NULL,
                     tolerance = 1e-10) {
    if (!is.numeric(g) || !length(g) || any(!is.finite(g)) ||
        !is.numeric(weights) ||
        length(weights) != length(g) || any(!is.finite(weights)) ||
        any(weights < 0) ||
        sum(weights) <= 0)
        stop("Invalid estimating values or weights.", call. = FALSE)
    keep <- weights > 0
    g <- g[keep]
    weights <- weights[keep]
    original.n <- sum(weights)
    if (!is.null(adjustment)) {
        .eq_scalar(adjustment, "adjustment", 0, Inf, TRUE)
        g <- c(g, -adjustment * sum(weights * g) / original.n)
        weights <- c(weights, 1)
    }
    .eq_scalar(tolerance, "tolerance", 0, .01, TRUE)
    if (all(g == 0))
        return(
            list(
                statistic = 0,
                lambda = 0,
                probabilities = weights / sum(weights),
                feasible = TRUE,
                residual = 0
            )
        )
    if (min(g) >= 0 || max(g) <= 0)
        return(
            list(
                statistic = Inf,
                lambda = NA_real_,
                probabilities = rep(NA_real_, length(g)),
                feasible = FALSE,
                residual = NA_real_
            )
        )
    scale <- max(abs(g))
    z <- g / scale
    score <- function(t)
        sum(weights * z / (1 + t * z))
    if (abs(sum(weights * z)) < 1e-14 * sum(weights))
        t <- 0
    else {
        left <- -1 / max(z)
        right <- -1 / min(z)
        # Root lies between zero and one pole, avoiding a huge irrelevant bracket.
        if (score(0) > 0) {
            lower <- 0
            upper <- right * (1 - 1e-13)
        } else {
            lower <- left * (1 - 1e-13)
            upper <- 0
        }
        if (!is.finite(score(lower)) || !is.finite(score(upper)) ||
            score(lower) < 0 || score(upper) > 0)
            stop("Could not bracket the EL dual root.", call. = FALSE)
        t <- stats::uniroot(score, c(lower, upper), tol = tolerance)$root
    }
    den <- 1 + t * z
    mass <- weights / (sum(weights) * den)
    list(
        statistic = max(0, 2 * sum(weights * log1p(t * z))),
        lambda = t / scale,
        probabilities = mass,
        feasible = TRUE,
        residual = sum(mass * g)
    )
}

.eq_ael_counts <- function(n, tau, adjustment, tolerance) {
    vapply(0:n, function(m)
        elq_dual(c(1 - tau, -tau), c(m, n - m), adjustment, tolerance)$statistic, numeric(1))
}

elq_bartlett <- function(x,
                         tau = .5,
                         bandwidth = NULL,
                         kernel = "biweight",
                         type = c("plugin", "limit")) {
    type <- match.arg(type)
    .eq_tau(tau)
    x <- .eq_sample(x)
    n <- length(x)
    if (type == "limit")
        return((1 - tau + tau^2) / (6 * tau * (1 - tau)))
    kernel <- .eq_kernels("chen_hall_bartlett", kernel)
    if (is.null(bandwidth))
        bandwidth <- elq_bandwidth(x, "chen_hall_bartlett")
    .eq_scalar(bandwidth, "bandwidth", 0, Inf, TRUE)
    q <- as.numeric(stats::quantile(x, tau, type = 1))
    g <- elq_kernel((q - x) / bandwidth, kernel, TRUE) - tau
    m2 <- mean(g^2)
    m3 <- mean(g^3)
    m4 <- mean(g^4)
    if (m2 <= .Machine$double.eps)
        stop("Degenerate moments for Bartlett correction.", call. = FALSE)
    m4 / (2 * m2^2) - m3^2 / (3 * m2^3)
}
