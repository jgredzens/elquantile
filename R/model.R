.eq_cdf_roots <- function(model, p) {
    # Exact piecewise-cubic root enumeration for the Zhou-Jing CDF analogue.
    x <- model$x
    h <- model$bandwidth
    knots <- sort(unique(c(x - h, x + h)))
    A <- (21 - 9 * sqrt(21)) / 8
    B <- (-3 + 3 * sqrt(21)) / 8
    roots <- numeric(0)
    for (i in seq_len(length(knots) - 1L)) {
        lo <- knots[i]
        hi <- knots[i + 1L]
        mid <- lo + (hi - lo) / 2
        radius <- (hi - lo) / 2
        z <- (mid - x) / h
        inside <- abs(z) < 1
        a <- z[inside]
        d <- radius / h
        n <- length(x)
        cf <- c((sum(z >= 1) + length(a) / 2 + B * sum(a) + A * sum(a^3) / 3) /
                    n - p,
                d * (B * length(a) + A * sum(a^2)) / n,
                A * d^2 * sum(a) / n,
                A * d^3 * length(a) / (3 * n)
        )
        while (length(cf) > 1L &&
               abs(tail(cf, 1L)) < 1e-15)
            cf <- head(cf, -1L)
        if (length(cf) == 1L) {
            if (abs(cf) < 1e-12)
                roots <- c(roots, lo, hi)
        } else {
            r <- polyroot(cf)
            good <- abs(Im(r)) < 1e-7 & Re(r) >= -1 - 1e-9 &
                Re(r) <= 1 + 1e-9
            roots <- c(roots, mid + radius * pmax(-1, pmin(1, Re(r[good]))))
        }
    }
    sort(unique(roots))
}

.eq_model <- function(x,
                      tau,
                      method,
                      bandwidth,
                      kernel,
                      adjustment,
                      bartlett,
                      control) {
    n <- length(x)
    smoothed <- method %in% c("chen_hall", "chen_hall_bartlett", "zhou_jing")
    if (!smoothed &&
        !is.null(bandwidth))
        stop("bandwidth is unused for this method.", call. = FALSE)
    if (smoothed) {
        if (is.null(bandwidth))
            bandwidth <- elq_bandwidth(x, method)
        .eq_scalar(bandwidth, "bandwidth", 0, Inf, TRUE)
    }
    if (method != "ael" &&
        !is.null(adjustment))
        stop("adjustment is only used for ael.", call. = FALSE)
    if (method == "ael" &&
        is.null(adjustment))
        adjustment <- max(1, log(n) / 2)
    if (!is.null(adjustment))
        .eq_scalar(adjustment, "adjustment", 0, Inf, TRUE)
    sample.q <- as.numeric(stats::quantile(x, tau, type = 1))
    model <- list(
        x = x,
        n = n,
        tau = tau,
        method = method,
        bandwidth = bandwidth,
        kernel = kernel,
        sample.quantile = sample.q,
        adjustment = adjustment
    )
    if (method %in% c("el", "ael")) {
        costs <- if (method == "el")
            .eq_binary(0:n, n, tau)
        else
            .eq_ael_counts(n, tau, adjustment, control$tolerance / 10)
        model$costs <- costs
        model$value <- function(theta)
            costs[findInterval(theta, x) + 1L]
        model$minimum <- min(costs)
        model$center <- sample.q
        model$knots <- x
        model$bounds <- if (method == "ael")
            c(-Inf, Inf)
        else
            range(x)
    } else if (method == "adimari") {
        probs <- (seq_len(n) - .5) / n
        cdf <- function(theta) {
            p <- stats::approx(x, probs, theta, rule = 2)$y
            p[theta < x[1L]] <- 0
            p[theta >= x[n]] <- 1
            p
        }
        model$cdf <- cdf
        model$value <- function(theta) {
            ans <- 2 * n * .eq_kl(cdf(theta), tau)
            ans[theta < x[1L] | theta >= x[n]] <- Inf
            ans
        }
        pbest <- max(probs[1L], min(probs[n], tau))
        model$minimum <- 2 * n * .eq_kl(pbest, tau)
        model$center <- as.numeric(stats::approx(probs, x, pbest, rule = 2)$y)
        model$knots <- x
        model$bounds <- range(x)
    } else {
        cdf <- function(theta)
            vapply(theta, function(q)
                mean(elq_kernel((q - x) / bandwidth, kernel, TRUE)), numeric(1))
        model$cdf <- cdf
        radius <- if (kernel == "gaussian")
            10
        else
            1
        model$bounds <- c(min(x) - radius * bandwidth, max(x) + radius * bandwidth)
        model$knots <- sort(unique(c(
            x - radius * bandwidth, x, x + radius * bandwidth
        )))
        if (method == "zhou_jing") {
            model$value <- function(theta) {
                p <- cdf(theta)
                ans <- 2 * n * .eq_kl(p, tau)
                ans[p <= 0 | p >= 1] <- Inf
                ans
            }
            roots <- .eq_cdf_roots(model, tau)
            if (!length(roots))
                stop("Failed to locate a Zhou-Jing CDF root.", call. = FALSE)
            model$center <- roots[which.min(abs(roots - sample.q))]
            model$roots <- roots
            model$knots <- sort(unique(c(model$knots, roots)))
            model$bartlett.factor <- NA_real_
        } else {
            factor <- if (method == "chen_hall_bartlett")
                elq_bartlett(x, tau, bandwidth, kernel, bartlett)
            else
                0
            multiplier <- 1 + factor / n
            if (!is.finite(multiplier) ||
                multiplier <= 0)
                stop("Invalid Bartlett scale.", call. = FALSE)
            model$value <- function(theta)
                vapply(theta, function(q)
                    elq_dual(
                        elq_kernel((q - x) / bandwidth, kernel, TRUE) - tau,
                        tolerance = control$tolerance / 10
                    )$statistic / multiplier, numeric(1))
            model$center <- stats::uniroot(function(q)
                cdf(q) - tau,
                model$bounds,
                tol = control$tolerance)$root
            model$bartlett.factor <- factor
        }
        model$minimum <- 0
    }
    model
}

.eq_models <- function(samples,
                       tau,
                       method,
                       bandwidth,
                       kernel,
                       adjustment,
                       bartlett,
                       control) {
    .eq_tau(tau)
    method <- .eq_method(method)
    kernel <- .eq_kernels(method, kernel)
    bartlett <- match.arg(bartlett, c("plugin", "limit"))
    if (!is.null(bandwidth) && (
        !is.numeric(bandwidth) ||
        !(length(bandwidth) %in% c(1L, length(samples))) ||
        any(!is.finite(bandwidth)) || any(bandwidth <= 0)
    ))
        stop("bandwidth must be positive, scalar or one per group.", call. = FALSE)
    if (!is.null(adjustment) && (
        !is.numeric(adjustment) ||
        !(length(adjustment) %in% c(1L, length(samples))) ||
        any(!is.finite(adjustment)) || any(adjustment <= 0)
    ))
        stop("adjustment must be positive, scalar or one per group.", call. = FALSE)
    h <- if (is.null(bandwidth))
        rep(list(NULL), length(samples))
    else
        as.list(rep(bandwidth, length.out = length(samples)))
    a <- if (is.null(adjustment))
        rep(list(NULL), length(samples))
    else
        as.list(rep(adjustment, length.out = length(samples)))
    Map(function(x, hh, aa)
        .eq_model(x, tau, method, hh, kernel, aa, bartlett, control),
        samples,
        h,
        a)
}

elq_statistic <- function(x,
                          theta,
                          tau = .5,
                          method = "el",
                          bandwidth = NULL,
                          kernel = NULL,
                          adjustment = NULL,
                          bartlett = "plugin",
                          control = elq_control()) {
    x <- .eq_sample(x)
    control <- .eq_control(control)
    if (!is.numeric(theta) ||
        !length(theta) || any(!is.finite(theta)))
        stop("theta must be a finite numeric vector.", call. = FALSE)
    model <- .eq_models(list(x),
                        tau,
                        method,
                        bandwidth,
                        kernel,
                        adjustment,
                        bartlett,
                        control)[[1L]]
    model$value(theta)
}

elq_profile <- function(x,
                        theta,
                        tau = .5,
                        method = "el",
                        bandwidth = NULL,
                        kernel = NULL,
                        adjustment = NULL,
                        bartlett = "plugin",
                        control = elq_control()) {
    structure(
        data.frame(
            theta = theta,
            statistic = elq_statistic(
                x,
                theta,
                tau,
                method,
                bandwidth,
                kernel,
                adjustment,
                bartlett,
                control
            )
        ),
        class = c("elq_profile", "data.frame"),
        tau = tau,
        method = method
    )
}
plot.elq_profile <- function(x, conf.level = .95, ...) {
    .eq_scalar(conf.level, "conf.level", 0, 1, TRUE)
    graphics::plot(
        x$theta,
        x$statistic,
        type = "l",
        xlab = "Quantile",
        ylab = "-2 log ratio",
        ...
    )
    graphics::abline(
        h = stats::qchisq(conf.level, 1),
        lty = 2,
        col = "red"
    )
    invisible(x)
}
