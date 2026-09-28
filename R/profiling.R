.eq_numerical_min <- function(fun, grid, tol, local = TRUE) {
    grid <- sort(unique(grid))
    values <- fun(grid)
    ok <- which(is.finite(values))
    if (!length(ok))
        return(list(value = Inf, at = NA_real_))
    i <- ok[which.min(values[ok])]
    ans <- list(value = values[i], at = grid[i])
    if (!local || length(grid) < 2L)
        return(ans)
    valleys <- ok[ok > 1L & ok < length(grid)]
    valleys <- valleys[values[valleys] <= values[valleys - 1L] &
                           values[valleys] <= values[valleys + 1L]]
    for (j in unique(c(valleys, i))) {
        interval <- grid[c(max(1L, j - 1L), min(length(grid), j + 1L))]
        if (diff(interval) <= 0)
            next
        fit <- suppressWarnings(stats::optimize(function(q) {
            w <- fun(q)
            if (is.finite(w))
                w
            else
                .Machine$double.xmax / 100
        }, interval, tol = tol))
        value <- fun(fit$minimum)
        if (is.finite(value) &&
            value < ans$value)
            ans <- list(value = value, at = fit$minimum)
    }
    ans
}

.eq_profile_groups <- function(models, offsets, control) {
    method <- models[[1L]]$method
    base <- sum(vapply(models, function(m)
        m$minimum, numeric(1)))
    fun <- function(theta) {
        z <- numeric(length(theta))
        for (j in seq_along(models))
            z <- z + models[[j]]$value(theta + offsets[j])
        z
    }
    knots <- sort(unique(unlist(
        Map(function(m, o)
            m$knots - o, models, offsets)
    )))
    if (length(knots) > control$max.cells)
        stop("Profile search exceeds max.cells.", call. = FALSE)
    if (method %in% c("el", "ael")) {
        # Include both exterior cells. Their actual thresholds do not affect counts.
        grid <- c(-Inf, knots)
        # Compare against shifted observations directly. Reconstructing theta+offset
        # can cross a sample jump by one floating-point unit and miss the best cell.
        values <- numeric(length(grid))
        for (j in seq_along(models))
            values <- values + models[[j]]$costs[findInterval(grid, models[[j]]$x -
                                                                  offsets[j]) + 1L]
        i <- which.min(values)
        return(
            list(
                statistic = max(0, values[i] - base),
                common = grid[i],
                baseline = base,
                constrained = values[i],
                optimizer = "Exact pooled CDF-cell enumeration",
                refinement.gap = 0
            )
        )
    }
    lower <- max(vapply(models, function(m)
        m$bounds[1L], numeric(1)) - offsets)
    upper <- min(vapply(models, function(m)
        m$bounds[2L], numeric(1)) - offsets)
    if (lower >= upper)
        return(
            list(
                statistic = Inf,
                common = NA_real_,
                baseline = base,
                constrained = Inf,
                optimizer = "Empty common support",
                refinement.gap = 0
            )
        )
    knots <- sort(unique(c(lower, knots[knots > lower &
                                            knots < upper], upper)))
    if (method == "adimari") {
        # On every open segment each interpolated CDF is affine and KL is convex.
        best <- list(value = Inf, at = NA_real_)
        for (i in seq_len(length(knots) - 1L)) {
            fit <- .eq_numerical_min(fun, c(knots[i], mean(knots[i:(i + 1L)]), knots[i +
                                                                                         1L]), control$tolerance)
            if (fit$value < best$value)
                best <- fit
        }
        gap <- 0
        optimizer <- "All interpolation segments; convex minimization within each"
    } else {
        centers <- vapply(models, function(m)
            m$center, numeric(1)) - offsets
        if (method != "zhou_jing") {
            lower <- max(lower, min(centers))
            upper <- min(upper, max(centers))
        }
        if (upper < lower)
            return(
                list(
                    statistic = Inf,
                    common = NA_real_,
                    baseline = base,
                    constrained = Inf,
                    optimizer = "No feasible center range",
                    refinement.gap = 0
                )
            )
        grid <- sort(unique(c(
            seq(lower, upper, length.out = control$grid.points),
            knots[knots >= lower &
                      knots <= upper],
            centers[centers >= lower & centers <= upper]
        )))
        best <- .eq_numerical_min(fun, grid, control$tolerance)
        gap <- NA_real_
        if (control$refine && length(grid) > 1L) {
            fine.grid <- sort(unique(c(grid, (
                head(grid, -1L) + tail(grid, -1L)
            ) / 2)))
            fine <- .eq_numerical_min(fun, fine.grid, control$tolerance)
            gap <- if (is.finite(fine$value) &&
                       is.finite(best$value))
                abs(best$value - fine$value)
            else
                if (identical(fine$value, best$value))
                    0
            else
                Inf
            if (fine$value < best$value)
                best <- fine
            if (gap > 1e-5 && control$warn)
                warning(
                    "Profile changed under grid refinement; inspect the numerical diagnostic.",
                    call. = FALSE
                )
        }
        optimizer <- "Knots, grid and local refinements; global optimum not certified"
    }
    list(
        statistic = max(0, best$value - base),
        common = best$at,
        baseline = base,
        constrained = best$value,
        optimizer = optimizer,
        refinement.gap = gap
    )
}
