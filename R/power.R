.eq_design <- function(n, quantiles, density, tau, contrast, null) {
    .eq_tau(tau)
    if (!is.numeric(quantiles) ||
        !length(quantiles) || any(!is.finite(quantiles)))
        stop("quantiles must be a finite numeric vector.", call. = FALSE)
    k <- length(quantiles)
    if (!is.numeric(n) ||
        !(length(n) %in% c(1L, k)) ||
        any(!is.finite(n)) || any(n < 2 | n != floor(n)))
        stop("n must contain integer sample sizes >= 2, scalar or one per group.",
             call. = FALSE)
    if (!is.numeric(density) || !(length(density) %in% c(1L, k)) ||
        any(!is.finite(density)) || any(density <= 0))
        stop("density must be positive, scalar or one density at each target quantile.",
             call. = FALSE)
    n <- rep(n, length.out = k)
    density <- rep(density, length.out = k)
    if (is.null(contrast))
        contrast <- if (k == 1L)
            matrix(1, 1, 1)
    else
        cbind(diag(k - 1L), -1)
    if (is.numeric(contrast) &&
        is.null(dim(contrast)))
        contrast <- matrix(contrast, nrow = 1L)
    if (!is.matrix(contrast) ||
        !is.numeric(contrast) || ncol(contrast) != k ||
        !nrow(contrast) ||
        any(!is.finite(contrast)) || qr(contrast)$rank != nrow(contrast))
        stop("contrast must be a finite, full-row-rank matrix with one column per group.",
             call. = FALSE)
    if (!is.numeric(null) ||
        !(length(null) %in% c(1L, nrow(contrast))) ||
        any(!is.finite(null)))
        stop("null must be scalar or one value per contrast.", call. = FALSE)
    null <- rep(null, length.out = nrow(contrast))
    effect <- drop(contrast %*% quantiles) - null
    variance <- tau * (1 - tau) / (n * density^2)
    covariance <- tcrossprod(sweep(contrast, 2L, sqrt(variance), "*"))
    list(
        n = n,
        quantiles = quantiles,
        density = density,
        tau = tau,
        contrast = contrast,
        null = null,
        effect = effect,
        covariance = covariance
    )
}

elq_power <- function(n,
                      quantiles,
                      density,
                      tau = .5,
                      contrast = NULL,
                      null = 0,
                      alpha = .05,
                      alternative = c("two.sided", "less", "greater")) {
    .eq_scalar(alpha, "alpha", 0, 1, TRUE)
    alternative <- match.arg(alternative)
    d <- .eq_design(n, quantiles, density, tau, contrast, null)
    df <- nrow(d$contrast)
    ncp <- max(0, drop(crossprod(
        d$effect, solve(d$covariance, d$effect)
    )))
    if (alternative == "two.sided") {
        power <- stats::pchisq(stats::qchisq(1 - alpha, df),
                               df,
                               ncp = ncp,
                               lower.tail = FALSE)
    } else {
        if (df != 1L)
            stop("Directional power requires a single contrast.", call. = FALSE)
        z <- d$effect / sqrt(d$covariance[1L, 1L])
        power <- if (alternative == "greater")
            stats::pnorm(z - stats::qnorm(1 - alpha))
        else
            stats::pnorm(stats::qnorm(alpha) - z)
    }
    structure(c(
        d,
        list(
            power = unname(power),
            alpha = alpha,
            df = df,
            ncp = ncp,
            alternative = alternative,
            method = "Local asymptotic quantile-contrast power",
            note = "Shared first-order approximation; does not distinguish finite-sample EL variants."
        )
    ), class = "elq_power")
}

elq_sample_size <- function(power = .8,
                            quantiles,
                            density,
                            ratio = rep(1, length(quantiles)),
                            tau = .5,
                            contrast = NULL,
                            null = 0,
                            alpha = .05,
                            alternative = c("two.sided", "less", "greater"),
                            min.n = 2L,
                            max.n = 1000000L) {
    .eq_scalar(power, "power", 0, 1, TRUE)
    alternative <- match.arg(alternative)
    min.n <- .eq_integer(min.n, "min.n", 2L)
    max.n <- .eq_integer(max.n, "max.n", min.n)
    k <- length(quantiles)
    if (!is.numeric(ratio) ||
        length(ratio) != k || any(!is.finite(ratio)) || any(ratio <= 0))
        stop("ratio must contain one positive allocation ratio per group.",
             call. = FALSE)
    ratio <- ratio / min(ratio)
    evaluate <- function(base)
        elq_power(ceiling(base * ratio),
                  quantiles,
                  density,
                  tau,
                  contrast,
                  null,
                  alpha,
                  alternative)
    first <- evaluate(min.n)
    if (first$power >= power)
        selected <- min.n
    else {
        high <- min.n
        while (high < max.n &&
               evaluate(high)$power < power)
            high <- min(max.n, 2 * high)
        if (evaluate(high)$power < power)
            stop(
                "Target power is not reached within max.n; check effect direction, densities, and search limit.",
                call. = FALSE
            )
        low <- min.n
        while (high - low > 1L) {
            mid <- floor(low + (high - low) / 2)
            if (evaluate(mid)$power >= power)
                high <- mid
            else
                low <- mid
        }
        selected <- high
    }
    ans <- evaluate(selected)
    structure(
        list(
            n = ans$n,
            total.n = sum(ans$n),
            base.n = selected,
            ratio = ratio,
            requested.power = power,
            power = ans$power,
            previous.power = if (selected > min.n)
                evaluate(selected - 1L)$power
            else
                NA_real_,
            alpha = alpha,
            tau = tau,
            planning = ans,
            method = "Smallest integer base size on the specified allocation path"
        ),
        class = "elq_sample_size"
    )
}

elq_density <- function(x, tau = .5, bandwidth = NULL) {
    x <- .eq_sample(x)
    .eq_tau(tau)
    if (is.null(bandwidth))
        bandwidth <- as.numeric(stats::IQR(x)) / (2 * stats::qnorm(.75)) * length(x)^(-.2)
    .eq_scalar(bandwidth, "bandwidth", 0, Inf, TRUE)
    q <- as.numeric(stats::quantile(x, tau, type = 1))
    value <- mean(stats::dnorm((q - x) / bandwidth)) / bandwidth
    structure(
        value,
        quantile = q,
        tau = tau,
        bandwidth = bandwidth,
        note = "Pilot KDE estimate; density uncertainty is not propagated into power."
    )
}

.eq_count_p <- function(n, tau, method, alternative, adjustment) {
    m <- 0:n
    if (method == "sign") {
        if (alternative == "greater")
            return(stats::pbinom(m, n, tau))
        if (alternative == "less")
            return(stats::pbinom(m - 1L, n, tau, lower.tail = FALSE))
        # Same probability ordering and relative tolerance as stats::binom.test,
        # evaluated for all counts in O(n log n) rather than n separate tests.
        mass <- stats::dbinom(m, n, tau)
        ordered <- sort(mass)
        indices <- findInterval(mass * (1 + 1e-7), ordered)
        return(pmin(1, c(0, cumsum(ordered))[indices + 1L]))
    }
    w <- if (method == "el")
        .eq_binary(m, n, tau)
    else
        .eq_ael_counts(n, tau, if (is.null(adjustment))
            max(1, log(n) / 2)
            else
                adjustment, 1e-11)
    if (alternative == "two.sided")
        return(stats::pchisq(w, 1, lower.tail = FALSE))
    direction <- sign(tau - m / n)
    z <- numeric(length(w))
    keep <- direction != 0
    z[keep] <- direction[keep] * sqrt(w[keep])
    if (alternative == "greater")
        stats::pnorm(z, lower.tail = FALSE)
    else
        stats::pnorm(z)
}

elq_power_exact <- function(n,
                            p,
                            tau = .5,
                            method = c("el", "ael", "sign"),
                            alpha = .05,
                            alternative = c("two.sided", "less", "greater"),
                            adjustment = NULL) {
    n <- .eq_integer(n, "n", 2L)
    .eq_tau(tau)
    .eq_scalar(p, "p", 0, 1)
    .eq_scalar(alpha, "alpha", 0, 1, TRUE)
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    if (n > 100000L)
        stop("Exact count enumeration is limited to n <= 100000.", call. = FALSE)
    if (!is.null(adjustment) &&
        method != "ael")
        stop("adjustment applies only to ael.", call. = FALSE)
    pv <- .eq_count_p(n, tau, method, alternative, adjustment)
    reject <- pv <= alpha
    structure(
        list(
            n = n,
            power = sum(stats::dbinom(0:n, n, p)[reject]),
            null.size = sum(stats::dbinom(0:n, n, tau)[reject]),
            p = p,
            tau = tau,
            alpha = alpha,
            alternative = alternative,
            method = paste("Exact binomial enumeration of", method),
            counts = data.frame(
                count = 0:n,
                p.value = pv,
                reject = reject
            ),
            note = "p is F_alternative(q0). Exact power enumeration does not make chi-square-calibrated EL an exact-level test."
        ),
        class = "elq_power"
    )
}

elq_sample_size_exact <- function(power = .8,
                                  p,
                                  tau = .5,
                                  method = c("el", "ael", "sign"),
                                  alpha = .05,
                                  alternative = c("two.sided", "less", "greater"),
                                  min.n = 2L,
                                  max.n = 2000L,
                                  adjustment = NULL,
                                  max.null.size = NULL) {
    .eq_scalar(power, "power", 0, 1, TRUE)
    method <- match.arg(method)
    alternative <- match.arg(alternative)
    min.n <- .eq_integer(min.n, "min.n", 2L)
    max.n <- .eq_integer(max.n, "max.n", min.n)
    if (!is.null(max.null.size))
        .eq_scalar(max.null.size, "max.null.size", 0, 1, TRUE)
    history <- list()
    selected <- NULL
    # Discrete power can oscillate: exhaustive ascending search is intentional.
    for (n in seq.int(min.n, max.n)) {
        z <- elq_power_exact(n, p, tau, method, alpha, alternative, adjustment)
        history[[length(history) + 1L]] <- data.frame(
            n = n,
            power = z$power,
            null.size = z$null.size
        )
        if (z$power >= power &&
            (is.null(max.null.size) || z$null.size <= max.null.size)) {
            selected <- z
            break
        }
    }
    if (is.null(selected))
        stop(
            "No sample size satisfies the power and null-size criteria within the search range.",
            call. = FALSE
        )
    structure(
        list(
            n = selected$n,
            total.n = selected$n,
            power = selected$power,
            requested.power = power,
            null.size = selected$null.size,
            alpha = alpha,
            tau = tau,
            history = do.call(rbind, history),
            planning = selected,
            method = "First qualifying integer n by exhaustive exact-count enumeration"
        ),
        class = "elq_sample_size"
    )
}

print.elq_power <- function(x, ...) {
    cat(x$method, "\n")
    cat("n:",
        paste(x$n, collapse = ", "),
        "\nPower:",
        format(x$power, digits = 5),
        "\n")
    if (!is.null(x$null.size))
        cat("Actual null rejection probability:",
            format(x$null.size, digits = 5),
            "\n")
    if (!is.null(x$mcse))
        cat("Monte Carlo SE:", format(x$mcse, digits = 4), "\n")
    if (!is.null(x$failures))
        cat("Failed replicates:", x$failures, "\n")
    cat(x$note, "\n")
    invisible(x)
}
print.elq_sample_size <- function(x, ...) {
    cat(x$method, "\n")
    cat(
        "Group sizes:",
        paste(x$n, collapse = ", "),
        "\nTotal:",
        x$total.n,
        "\nPower:",
        format(x$power, digits = 5),
        "\n"
    )
    invisible(x)
}
