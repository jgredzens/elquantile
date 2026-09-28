.eq_wilson <- function(r, n, level = .95) {
    if (n == 0)
        return(c(lower = NA_real_, upper = NA_real_))
    z <- stats::qnorm(1 - (1 - level) / 2)
    p <- r / n
    center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    radius <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(lower = center - radius, upper = center + radius)
}

elq_power_sim <- function(n,
                          generators,
                          tau = .5,
                          null = 0,
                          paired = FALSE,
                          method = "el",
                          nsim = 1000L,
                          alpha = .05,
                          alternative = c("two.sided", "less", "greater"),
                          seed = NULL,
                          keep = TRUE,
                          ...) {
    .eq_flag(paired, "paired")
    .eq_flag(keep, "keep")
    .eq_tau(tau)
    nsim <- .eq_integer(nsim, "nsim", 2L)
    .eq_scalar(alpha, "alpha", 0, 1, TRUE)
    alternative <- match.arg(alternative)
    method <- .eq_method(method)
    dots <- list(...)
    if (any(
        names(dots) %in% c(
            "x",
            "y",
            "group",
            "tau",
            "null",
            "paired",
            "method",
            "alternative",
            "conf.int"
        )
    ))
    stop("Do not override simulation data or test arguments through ...",
         call. = FALSE)
    if (paired) {
        n <- .eq_integer(n, "n", 2L)
        if (!is.function(generators))
            stop("Paired generators must be a function(n) returning an n by 2 matrix.",
                 call. = FALSE)
    } else {
        if (is.function(generators))
            generators <- list(generators)
        if (!is.list(generators) ||
            !length(generators) ||
            !all(vapply(generators, is.function, logical(1))))
            stop(
                "Use a generator function for one sample or a list of functions for independent groups.",
                call. = FALSE
            )
        k <- length(generators)
        if (!is.numeric(n) ||
            !(length(n) %in% c(1L, k)) ||
            any(!is.finite(n)) || any(n < 2 | n != floor(n)))
            stop("Invalid group sizes.", call. = FALSE)
        n <- rep(n, length.out = k)
    }
    run <- function() {
        rows <- vector("list", nsim)
        for (i in seq_len(nsim)) {
            warnings <- character(0)
            ans <- tryCatch(
                withCallingHandlers({
                    if (paired) {
                        d <- generators(n)
                        if (!is.matrix(d) ||
                            !is.numeric(d) || nrow(d) != n || ncol(d) != 2L)
                            stop("Paired generator returned an invalid matrix.")
                        args <- list(x = d[, 1L],
                                     y = d[, 2L],
                                     paired = TRUE)
                    } else {
                        samples <- Map(function(fun, nn) {
                            z <- fun(nn)
                            if (!is.numeric(z) ||
                                length(z) != nn ||
                                !is.null(dim(z)))
                                stop("Generator returned the wrong shape.")
                            z
                        }, generators, n)
                        args <- list(x = if (length(samples) == 1L)
                            samples[[1L]]
                            else
                                samples)
                    }
                    fit <- do.call(elq_test, c(
                        args,
                        list(
                            tau = tau,
                            null = null,
                            method = method,
                            alternative = alternative,
                            conf.int = FALSE
                        ),
                        dots
                    ))
                    if (!is.finite(fit$p.value) ||
                        fit$p.value < 0 || fit$p.value > 1)
                        stop("Invalid p-value.")
                    fit
                }, warning = function(w) {
                    warnings <<- c(warnings, conditionMessage(w))
                    invokeRestart("muffleWarning")
                }),
                error = function(e)
                    e
            )
            ok <- !inherits(ans, "error")
            rows[[i]] <- data.frame(
                replicate = i,
                p.value = if (ok)
                    ans$p.value
                else
                    NA_real_,
                statistic = if (ok)
                    unname(ans$statistic)
                else
                    NA_real_,
                reject = if (ok)
                    ans$p.value <= alpha
                else
                    NA,
                status = if (ok)
                    "ok"
                else
                    "failed",
                error = if (ok)
                    ""
                else
                    conditionMessage(ans),
                warnings = paste(unique(warnings), collapse = " | "),
                refinement.gap = if (ok)
                    ans$refinement.gap
                else
                    NA_real_,
                stringsAsFactors = FALSE
            )
        }
        do.call(rbind, rows)
    }
    raw <- .eq_seed(seed, run)
    good <- raw$status == "ok"
    valid <- sum(good)
    rejected <- sum(raw$reject[good])
    power <- if (valid)
        rejected / valid
    else
        NA_real_
    result <- list(
        n = n,
        nsim = nsim,
        valid = valid,
        failures = nsim - valid,
        power = power,
        mcse = if (valid)
            sqrt(power * (1 - power) / valid)
        else
            NA_real_,
        mc.interval = .eq_wilson(rejected, valid),
        failure.bounds = c(
            lower = rejected / nsim,
            upper = (rejected + nsim - valid) / nsim
        ),
        tau = tau,
        alpha = alpha,
        paired = paired,
        method.id = method,
        method = paste("Monte Carlo rejection probability:", method),
        note = "Power is conditional on successful replicates; inspect failures and worst-case bounds. It estimates size when the generator satisfies the null.",
        seed = seed
    )
    if (keep)
        result$replicates <- raw
    structure(result, class = "elq_power")
}

elq_sample_size_sim <- function(n.grid,
                                generators,
                                power = .8,
                                ratio = NULL,
                                paired = FALSE,
                                criterion = c("estimate", "lower"),
                                nsim = 1000L,
                                seed = NULL,
                                ...) {
    criterion <- match.arg(criterion)
    .eq_scalar(power, "power", 0, 1, TRUE)
    if (!is.numeric(n.grid) ||
        !length(n.grid) || any(!is.finite(n.grid)) ||
        any(n.grid < 2 |
            n.grid != floor(n.grid)) || anyDuplicated(n.grid))
        stop("n.grid must contain distinct integer base sizes >= 2.", call. = FALSE)
    n.grid <- sort(n.grid)
    k <- if (paired ||
             is.function(generators))
        1L
    else
        length(generators)
    if (is.null(ratio))
        ratio <- rep(1, k)
    if (!is.numeric(ratio) ||
        length(ratio) != k || any(!is.finite(ratio)) || any(ratio <= 0))
        stop("Invalid allocation ratio.", call. = FALSE)
    ratio <- ratio / min(ratio)
    run <- function()
        lapply(n.grid, function(b)
            elq_power_sim(
                ceiling(b * ratio),
                generators,
                paired = paired,
                nsim = nsim,
                seed = NULL,
                ...
            ))
    fits <- .eq_seed(seed, run)
    table <- do.call(rbind, lapply(seq_along(fits), function(i) {
        z <- fits[[i]]
        data.frame(
            base.n = n.grid[i],
            total.n = sum(z$n),
            power = z$power,
            mcse = z$mcse,
            mc.lower = z$mc.interval[1L],
            mc.upper = z$mc.interval[2L],
            failures = z$failures,
            row.names = NULL
        )
    }))
    score <- if (criterion == "estimate")
        table$power
    else
        table$mc.lower
    hits <- which(!is.na(score) &
                      score >= power & table$failures == 0)
    selected <- if (length(hits))
        hits[1L]
    else
        NA_integer_
    list(
        selected.n = if (is.na(selected))
            NULL
        else
            fits[[selected]]$n,
        table = table,
        simulations = fits,
        criterion = criterion,
        requested.power = power,
        note = "Selection is restricted to this grid. Pointwise Monte Carlo intervals do not give simultaneous or post-selection guarantees; confirm the selected design independently."
    )
}
