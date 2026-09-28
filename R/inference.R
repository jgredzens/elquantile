.eq_pvalue <- function(w, df, direction, alternative) {
    if (alternative == "two.sided")
        return(stats::pchisq(w, df, lower.tail = FALSE))
    if (df != 1L)
        stop("One-sided alternatives require one sample or two groups.",
             call. = FALSE)
    z <- if (direction == 0)
        0
    else
        sign(direction) * sqrt(w)
    if (alternative == "greater")
        stats::pnorm(z, lower.tail = FALSE)
    else
        stats::pnorm(z)
}

elq_test <- function(x,
                     y = NULL,
                     group = NULL,
                     tau = .5,
                     null = 0,
                     paired = FALSE,
                     method = "el",
                     alternative = c("two.sided", "less", "greater"),
                     conf.int = FALSE,
                     conf.level = .95,
                     bandwidth = NULL,
                     kernel = NULL,
                     adjustment = NULL,
                     bartlett = "plugin",
                     control = elq_control()) {
    method <- .eq_method(method)
    alternative <- match.arg(alternative)
    control <- .eq_control(control)
    .eq_flag(paired, "paired")
    .eq_flag(conf.int, "conf.int")
    .eq_scalar(conf.level, "conf.level", 0, 1, TRUE)
    .eq_scalar(null, "null")
    paired.target <- FALSE
    if (paired) {
        if (!is.null(group) ||
            is.list(x) || is.null(y) || length(x) != length(y))
            stop("Paired inference needs equal-length numeric x and y.",
                 call. = FALSE)
        if (!is.numeric(x) ||
            !is.numeric(y) || any(!is.finite(x)) || any(!is.finite(y)))
            stop("Paired observations must be finite numeric values.",
                 call. = FALSE)
        x <- x - y
        y <- NULL
        paired.target <- TRUE
    }
    if (!is.null(group) || (is.list(x) && !is.data.frame(x))) {
        if (!is.null(y))
            stop("Do not combine grouped input with y.", call. = FALSE)
        samples <- .eq_groups(x, group)
    } else if (is.null(y))
        samples <- list(x = .eq_sample(x))
    else
        samples <- list(x = .eq_sample(x), y = .eq_sample(y, "y"))
    k <- length(samples)
    if (k > 2L && (null != 0 || alternative != "two.sided"))
        stop(
            "k-sample inference tests equality; use null=0 and alternative='two.sided'.",
            call. = FALSE
        )
    if (k > 2L &&
        conf.int)
        stop("Use elq_posthoc() for group confidence intervals.", call. = FALSE)
    models <- .eq_models(samples,
                         tau,
                         method,
                         bandwidth,
                         kernel,
                         adjustment,
                         bartlett,
                         control)
    .eq_notice(samples, tau, method, control)
    estimates <- vapply(models, function(m)
        m$sample.quantile, numeric(1))
    if (k == 1L) {
        fit <- list(
            statistic = models[[1L]]$value(null),
            common = NA_real_,
            baseline = 0,
            constrained = NA_real_,
            optimizer = "Pointwise evaluation",
            refinement.gap = 0
        )
        df <- 1L
        direction <- estimates[1L] - null
        target <- if (paired.target)
            "quantile of the paired differences x - y"
        else
            "population quantile"
        estimate <- c(quantile = unname(estimates[1L]))
    } else {
        offsets <- c(if (k == 2L)
            null
            else
                0, rep(0, k - 1L))
        fit <- .eq_profile_groups(models, offsets, control)
        df <- k - 1L
        direction <- if (k == 2L)
            estimates[1L] - estimates[2L] - null
        else
            NA_real_
        target <- if (k == 2L)
            "difference of independent population quantiles q_x - q_y"
        else
            "equality of independent population quantiles"
        estimate <- if (k == 2L)
            stats::setNames(unname(estimates[1L] - estimates[2L]), "quantile difference")
        else
            estimates
    }
    result <- list(
        statistic = c(W = fit$statistic),
        parameter = c(df = df),
        p.value = .eq_pvalue(fit$statistic, df, direction, alternative),
        alternative = alternative,
        null.value = if (k <= 2L)
            c(quantile = null)
        else
            NULL,
        estimate = estimate,
        method = paste(
            elq_methods()$description[match(method, elq_methods()$method)],
            if (k == 1L)
                "(one sample)"
            else
                paste0("(", k, " independent groups)")
        ),
        data.name = paste(names(samples), collapse = ", "),
        target = target,
        method.id = method,
        tau = tau,
        n = lengths(samples),
        paired = paired.target,
        calibration = "First-order chi-square; one-sided signed-root normal when requested",
        bandwidth = vapply(models, function(m)
            if (is.null(m$bandwidth))
                NA_real_
            else
                m$bandwidth, numeric(1)),
        bartlett.factor = vapply(models, function(m)
            if (is.null(m$bartlett.factor))
                NA_real_
            else
                m$bartlett.factor, numeric(1)),
        group.estimates = estimates,
        common.quantile = fit$common,
        baseline = fit$baseline,
        optimizer = fit$optimizer,
        refinement.gap = fit$refinement.gap
    )
    if (conf.int) {
        result$conf.int <- if (k == 1L)
            .eq_one_ci(models[[1L]], conf.level, control)
        else
            .eq_difference_ci(models, conf.level, "bonferroni", control)
        result$interval.note <- "Two-sided interval; independent differences use Bonferroni marginal hulls."
    }
    structure(result, class = c("elq_test", "htest"))
}

elq_one_sample <- function(x, q0 = 0, ...)
    elq_test(x, null = q0, ...)
elq_two_sample <- function(x,
                           y,
                           delta = 0,
                           paired = FALSE,
                           ...)
    elq_test(x, y, null = delta, paired = paired, ...)
elq_k_sample <- function(x, group = NULL, ...) {
    samples <- .eq_groups(x, group)
    elq_test(samples, ...)
}
elq_anova <- function(formula, data, ...) {
    if (!inherits(formula, "formula") ||
        length(formula) != 3L || !is.symbol(formula[[3L]]))
        stop("Use response ~ group with one grouping variable.", call. = FALSE)
    mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
    if (ncol(mf) != 2L)
        stop("Use one response and one group variable.", call. = FALSE)
    elq_k_sample(mf[[1L]], mf[[2L]], ...)
}
print.elq_test <- function(x, ...) {
    print(structure(x, class = "htest"), ...)
    cat("Target:", x$target, "\n")
    if (!is.null(x$interval.note))
        cat(x$interval.note, "\n")
    invisible(x)
}

elq_posthoc <- function(x,
                        group = NULL,
                        tau = .5,
                        method = "el",
                        p.adjust.method = c("holm", "bonferroni", "none"),
                        conf.int = FALSE,
                        conf.level = .95,
                        bandwidth = NULL,
                        kernel = NULL,
                        adjustment = NULL,
                        bartlett = "plugin",
                        control = elq_control()) {
    p.adjust.method <- match.arg(p.adjust.method)
    control <- .eq_control(control)
    .eq_flag(conf.int, "conf.int")
    .eq_scalar(conf.level, "conf.level", 0, 1, TRUE)
    if (!is.numeric(tau) ||
        !length(tau) ||
        any(!is.finite(tau)) ||
        any(tau <= 0 | tau >= 1) || anyDuplicated(tau))
        stop("tau must contain distinct probabilities in (0,1).", call. = FALSE)
    samples <- .eq_groups(x, group)
    pairs <- utils::combn(seq_along(samples), 2L)
    M <- ncol(pairs) * length(tau)
    rows <- list()
    level <- 1 - (1 - conf.level) / (2 * M)
    if (level >= 1)
        stop("Requested simultaneous interval level rounds to one.",
             call. = FALSE)
    for (tt in tau) {
        models <- .eq_models(samples,
                             tt,
                             method,
                             bandwidth,
                             kernel,
                             adjustment,
                             bartlett,
                             control)
        ci <- if (conf.int)
            lapply(models,
                   .eq_one_ci,
                   level = level,
                   control = control)
        else
            NULL
        for (i in seq_len(ncol(pairs))) {
            ab <- pairs[, i]
            fit <- .eq_profile_groups(models[ab], c(0, 0), control)
            row <- data.frame(
                group1 = names(samples)[ab[1L]],
                group2 = names(samples)[ab[2L]],
                tau = tt,
                estimate = models[[ab[1L]]]$sample.quantile - models[[ab[2L]]]$sample.quantile,
                statistic = fit$statistic,
                p.value = stats::pchisq(fit$statistic, 1, lower.tail = FALSE),
                refinement.gap = fit$refinement.gap,
                stringsAsFactors = FALSE
            )
            if (conf.int) {
                row$conf.low <- unname(ci[[ab[1L]]][1L] - ci[[ab[2L]]][2L])
                row$conf.high <- unname(ci[[ab[1L]]][2L] - ci[[ab[2L]]][1L])
            }
            rows[[length(rows) + 1L]] <- row
        }
    }
    out <- do.call(rbind, rows)
    out$p.adj <- stats::p.adjust(out$p.value, p.adjust.method)
    attr(out, "family.size") <- M
    attr(out, "p.adjust.method") <- p.adjust.method
    attr(out, "calibration") <- "Fixed-family first-order inference; all pairs and tau values form one family"
    if (conf.int)
        attr(out, "interval.method") <- "Bonferroni differences of marginal confidence hulls; not Holm inversion"
    out
}

elq_sign_test <- function(x,
                          q0 = 0,
                          tau = .5,
                          alternative = c("two.sided", "less", "greater"),
                          conf.level = .95) {
    x <- .eq_sample(x)
    .eq_scalar(q0, "q0")
    .eq_tau(tau)
    alternative <- match.arg(alternative)
    .eq_scalar(conf.level, "conf.level", 0, 1, TRUE)
    if (any(x == q0))
        stop(
            "An observation equals q0; continuous-null sign inference is ambiguous for this dataset.",
            call. = FALSE
        )
    m <- sum(x <= q0)
    n <- length(x)
    count.alt <- switch(
        alternative,
        less = "greater",
        greater = "less",
        two.sided = "two.sided"
    )
    z <- stats::binom.test(m, n, tau, alternative = count.alt)
    tail <- (1 - conf.level) / 2
    lower.index <- as.integer(stats::qbinom(tail, n, tau))
    upper.index <- as.integer(stats::qbinom(1 - tail, n, tau)) + 1L
    ci <- c(lower = if (lower.index < 1L)
        - Inf
        else
            x[lower.index], upper = if (upper.index > n)
                Inf
        else
            x[upper.index])
    attr(ci, "conf.level") <- conf.level
    structure(
        list(
            statistic = c(S = m),
            parameter = c(n = n),
            p.value = z$p.value,
            estimate = c(quantile = as.numeric(stats::quantile(x, tau, type = 1))),
            null.value = c(quantile = q0),
            alternative = alternative,
            method = "Exact binomial quantile sign test",
            data.name = "x",
            conf.int = ci,
            interval.note = "Two-sided equal-tail order-statistic interval",
            target = "population quantile",
            tau = tau
        ),
        class = c("elq_test", "htest")
    )
}
