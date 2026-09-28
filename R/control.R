# Input validation and numerical settings are separate.
.eq_scalar <- function(x,
                       name,
                       lower = -Inf,
                       upper = Inf,
                       open = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
    if (ok)
        ok <- if (open)
            x > lower && x < upper
    else
        x >= lower && x <= upper
    if (!ok)
        stop(name,
             " must be one finite number in the specified range.",
             call. = FALSE)
    invisible(x)
}
.eq_integer <- function(x, name, lower = 1L) {
    .eq_scalar(x, name, lower, .Machine$integer.max)
    if (x != floor(x))
        stop(name, " must be an integer.", call. = FALSE)
    as.integer(x)
}
.eq_tau <- function(tau)
    .eq_scalar(tau, "tau", 0, 1, TRUE)
.eq_flag <- function(x, name) {
    if (!is.logical(x) ||
        length(x) != 1L ||
        is.na(x))
        stop(name, " must be TRUE or FALSE.", call. = FALSE)
}
.eq_sample <- function(x, name = "x") {
    if (!is.numeric(x) ||
        !is.null(dim(x)) || length(x) < 2L || any(!is.finite(x)))
        stop(name,
             " must contain at least two finite numeric observations.",
             call. = FALSE)
    if (anyDuplicated(x))
        stop(
            name,
            " contains ties; the continuous-population methods require distinct observations.",
            call. = FALSE
        )
    sort(x)
}
.eq_groups <- function(x, group = NULL) {
    if (is.list(x) && !is.data.frame(x)) {
        if (!is.null(group))
            stop("Do not supply group with a list.", call. = FALSE)
        z <- x
    } else {
        if (is.null(group) || length(x) != length(group) || anyNA(group))
            stop("Supply a list of samples, or x and complete group labels.",
                 call. = FALSE)
        z <- split(x, factor(group), drop = TRUE)
    }
    if (length(z) < 2L)
        stop("At least two groups are required.", call. = FALSE)
    if (is.null(names(z)))
        names(z) <- paste0("group", seq_along(z))
    if (anyNA(names(z)) ||
        any(!nzchar(names(z))) || anyDuplicated(names(z)))
        stop("Group names must be nonempty and unique.", call. = FALSE)
    lapply(z, .eq_sample)
}

elq_control <- function(grid.points = 81L,
                        refine = TRUE,
                        tolerance = 1e-9,
                        max.cells = 2e6,
                        warn = TRUE) {
    .eq_integer(grid.points, "grid.points", 21L)
    .eq_flag(refine, "refine")
    .eq_flag(warn, "warn")
    .eq_scalar(tolerance, "tolerance", 0, 0.01, TRUE)
    .eq_scalar(max.cells, "max.cells", 1)
    structure(
        list(
            grid.points = as.integer(grid.points),
            refine = refine,
            tolerance = tolerance,
            max.cells = max.cells,
            warn = warn
        ),
        class = "elq_control"
    )
}
.eq_control <- function(control) {
    if (!inherits(control, "elq_control"))
        stop("Use elq_control() to construct control.", call. = FALSE)
    control
}
.eq_notice <- function(samples, tau, method, control) {
    if (control$warn && min(lengths(samples)) * min(tau, 1 - tau) < 5)
        warning(
            "Sparse quantile tail: asymptotic calibration may be inaccurate; consider exact or simulation-based assessment.",
            call. = FALSE
        )
    invisible(NULL)
}
.eq_seed <- function(seed, fun) {
    if (is.null(seed))
        return(fun())
    .eq_integer(seed, "seed", 0L)
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (existed)
        old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        if (existed)
            assign(".Random.seed", old, envir = .GlobalEnv)
        else
            if (exists(".Random.seed",
                       envir = .GlobalEnv,
                       inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
    fun()
}
.eq_interval <- function(lower, upper, level, method, components = NULL) {
    z <- c(lower = lower, upper = upper)
    attr(z, "conf.level") <- level
    attr(z, "method") <- method
    attr(z, "components") <- components
    class(z) <- c("elq_interval", "numeric")
    z
}
