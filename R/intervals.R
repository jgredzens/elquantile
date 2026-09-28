.eq_merge <- function(lower, upper) {
    if (!length(lower))
        return(data.frame(lower = numeric(0), upper = numeric(0)))
    o <- order(lower, upper)
    lower <- lower[o]
    upper <- upper[o]
    rows <- list()
    lo <- lower[1L]
    hi <- upper[1L]
    if (length(lower) > 1L)
        for (i in 2:length(lower)) {
            if (lower[i] <= hi)
                hi <- max(hi, upper[i])
            else {
                rows[[length(rows) + 1L]] <- c(lo, hi)
                lo <- lower[i]
                hi <- upper[i]
            }
        }
    rows[[length(rows) + 1L]] <- c(lo, hi)
    z <- as.data.frame(do.call(rbind, rows))
    names(z) <- c("lower", "upper")
    z
}
.eq_probability_limits <- function(n, tau, cutoff) {
    objective <- function(p)
        2 * n * .eq_kl(p, tau) - cutoff
    lower <- if (objective(0) <= 0)
        0
    else
        stats::uniroot(objective, c(0, tau), tol = 1e-12)$root
    upper <- if (objective(1) <= 0)
        1
    else
        stats::uniroot(objective, c(tau, 1), tol = 1e-12)$root
    c(lower, upper)
}
.eq_boundary <- function(fun, outside, inside, cutoff, tolerance) {
    for (i in seq_len(100L)) {
        mid <- outside + (inside - outside) / 2
        if (fun(mid) <= cutoff)
            inside <- mid
        else
            outside <- mid
        if (abs(inside - outside) <= tolerance * max(1, abs(mid)))
            break
    }
    inside
}
.eq_one_ci <- function(model, level, control) {
    cutoff <- stats::qchisq(level, 1)
    x <- model$x
    method <- model$method
    if (method %in% c("el", "ael")) {
        edges <- c(-Inf, x, Inf)
        keep <- model$costs <= cutoff
        components <- .eq_merge(head(edges, -1L)[keep], tail(edges, -1L)[keep])
    } else if (method == "adimari") {
        p <- .eq_probability_limits(model$n, model$tau, cutoff)
        probs <- (seq_len(model$n) - .5) / model$n
        p <- c(max(p[1L], probs[1L]), min(p[2L], tail(probs, 1L)))
        if (p[1L] > p[2L])
            components <- .eq_merge(numeric(0), numeric(0))
        else {
            endpoints <- stats::approx(probs, x, p, rule = 2)$y
            components <- .eq_merge(endpoints[1L], endpoints[2L])
        }
    } else if (method == "zhou_jing") {
        p <- .eq_probability_limits(model$n, model$tau, cutoff)
        breaks <- sort(unique(
            c(
                model$bounds,
                model$knots,
                .eq_cdf_roots(model, p[1L]),
                .eq_cdf_roots(model, p[2L]),
                .eq_cdf_roots(model, 0),
                .eq_cdf_roots(model, 1)
            )
        ))
        lo <- head(breaks, -1L)
        hi <- tail(breaks, -1L)
        keep <- model$value(lo + (hi - lo) / 2) <= cutoff
        components <- .eq_merge(lo[keep], hi[keep])
    } else {
        lower <- .eq_boundary(model$value,
                              model$bounds[1L],
                              model$center,
                              cutoff,
                              control$tolerance)
        upper <- .eq_boundary(model$value,
                              model$bounds[2L],
                              model$center,
                              cutoff,
                              control$tolerance)
        components <- .eq_merge(lower, upper)
    }
    if (!nrow(components))
        return(.eq_interval(
            NA_real_,
            NA_real_,
            level,
            "Empty confidence set",
            components
        ))
    .eq_interval(
        min(components$lower),
        max(components$upper),
        level,
        "Closed hull of chi-square-calibrated confidence set; components retained",
        components
    )
}
.eq_difference_ci <- function(models, level, interval, control) {
    if (interval == "bonferroni") {
        a <- .eq_one_ci(models[[1L]], 1 - (1 - level) / 2, control)
        b <- .eq_one_ci(models[[2L]], 1 - (1 - level) / 2, control)
        return(
            .eq_interval(
                unname(a[1L] - b[2L]),
                unname(a[2L] - b[1L]),
                level,
                "Bonferroni difference of two marginal confidence hulls; not inversion of the difference test"
            )
        )
    }
    if (!(models[[1L]]$method %in% c("el", "ael")))
        stop(
            "Profile difference intervals currently support el and ael; use interval='bonferroni' for all methods.",
            call. = FALSE
        )
    a <- models[[1L]]
    b <- models[[2L]]
    if ((a$n + 1) * (b$n + 1) > control$max.cells)
        stop("Difference CI exceeds max.cells.", call. = FALSE)
    ea <- c(-Inf, a$x, Inf)
    eb <- c(-Inf, b$x, Inf)
    cutoff <- stats::qchisq(level, 1) + a$minimum + b$minimum
    lower <- Inf
    upper <- -Inf
    for (i in seq_along(a$costs)) {
        keep <- which(a$costs[i] + b$costs <= cutoff)
        if (length(keep)) {
            lower <- min(lower, ea[i] - eb[keep + 1L])
            upper <- max(upper, ea[i + 1L] - eb[keep])
        }
    }
    if (lower > upper)
        lower <- upper <- NA_real_
    .eq_interval(lower,
                 upper,
                 level,
                 "Closed hull of exact step-profile difference confidence set")
}

elq_confint <- function(x,
                        y = NULL,
                        tau = .5,
                        paired = FALSE,
                        method = "el",
                        conf.level = .95,
                        interval = c("bonferroni", "profile"),
                        bandwidth = NULL,
                        kernel = NULL,
                        adjustment = NULL,
                        bartlett = "plugin",
                        control = elq_control()) {
    .eq_flag(paired, "paired")
    .eq_scalar(conf.level, "conf.level", 0, 1, TRUE)
    interval <- match.arg(interval)
    control <- .eq_control(control)
    if (paired) {
        if (is.null(y) ||
            length(x) != length(y))
            stop("Paired x and y must have equal lengths.", call. = FALSE)
        if (!is.numeric(x) ||
            !is.numeric(y) || any(!is.finite(x)) || any(!is.finite(y)))
            stop("Paired observations must be finite numeric values.",
                 call. = FALSE)
        x <- x - y
        y <- NULL
    }
    samples <- if (is.null(y))
        list(.eq_sample(x))
    else
        list(.eq_sample(x), .eq_sample(y, "y"))
    models <- .eq_models(samples,
                         tau,
                         method,
                         bandwidth,
                         kernel,
                         adjustment,
                         bartlett,
                         control)
    if (length(models) == 1L)
        .eq_one_ci(models[[1L]], conf.level, control)
    else
        .eq_difference_ci(models, conf.level, interval, control)
}
