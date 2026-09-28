# Method metadata here, implementations in model.R.
elq_methods <- function() {
    data.frame(
        method = c(
            "el",
            "adimari",
            "chen_hall",
            "chen_hall_bartlett",
            "zhou_jing",
            "ael"
        ),
        description = c(
            "Unsmoothed binary empirical likelihood",
            "Adimari midpoint interpolation in Bernoulli KL",
            "Chen-Hall smoothed estimating equations",
            "Chen-Hall with plug-in or limiting Bartlett factor",
            "Zhou-Jing signed-kernel substitution in Bernoulli KL",
            "Chen-Variyath-Abraham pseudo-observation adjustment"
        ),
        one.sample = TRUE,
        paired = TRUE,
        independent.two.sample = TRUE,
        k.sample = TRUE,
        group.profiling = c(
            "exact cells",
            "piecewise convex",
            "numerical",
            "numerical",
            "numerical, possibly nonconvex",
            "exact cells"
        ),
        reference = c(
            "Owen (1988, 1990)",
            "Adimari (1998)",
            "Chen and Hall (1993)",
            "Chen and Hall (1993), Section 4",
            "Zhou and Jing (2003)",
            "Chen, Variyath and Abraham (2008)"
        ),
        stringsAsFactors = FALSE
    )
}
.eq_method <- function(method)
    match.arg(method, elq_methods()$method)

elq_bandwidth <- function(x, method = "chen_hall", scale = NULL) {
    method <- .eq_method(method)
    x <- .eq_sample(x)
    if (!(method %in% c("chen_hall", "chen_hall_bartlett", "zhou_jing")))
        stop("This method does not use a bandwidth.", call. = FALSE)
    if (is.null(scale))
        scale <- as.numeric(stats::IQR(x, type = 7)) / (2 * stats::qnorm(.75))
    .eq_scalar(scale, "scale", 0, Inf, TRUE)
    exponent <- if (method == "chen_hall_bartlett")
        .75
    else
        .5
    scale * length(x)^(-exponent)
}

elq_kernel <- function(u,
                       kernel = c("biweight", "epanechnikov", "gaussian", "zhou_jing"),
                       cumulative = FALSE) {
    kernel <- match.arg(kernel)
    .eq_flag(cumulative, "cumulative")
    if (!is.numeric(u) ||
        anyNA(u))
        stop("u must be numeric without missing values.", call. = FALSE)
    if (kernel == "gaussian")
        return(if (cumulative)
            stats::pnorm(u)
            else
                stats::dnorm(u))
    z <- pmax(-1, pmin(1, u))
    if (kernel == "biweight") {
        ans <- if (cumulative)
            .5 + 15 / 16 * (z - 2 * z^3 / 3 + z^5 / 5)
        else
            15 / 16 * (1 - z^2)^2
    } else if (kernel == "epanechnikov") {
        ans <- if (cumulative)
            .5 + .75 * (z - z^3 / 3)
        else
            .75 * (1 - z^2)
    } else {
        A <- (21 - 9 * sqrt(21)) / 8
        B <- (-3 + 3 * sqrt(21)) / 8
        ans <- if (cumulative)
            .5 + B * z + A * z^3 / 3
        else
            A * z^2 + B
    }
    if (cumulative) {
        ans[u <= -1] <- 0
        ans[u >= 1] <- 1
    } else
        ans[abs(u) > 1] <- 0
    ans
}
.eq_kernels <- function(method, kernel) {
    if (method == "zhou_jing") {
        if (!is.null(kernel) && kernel != "zhou_jing")
            stop("zhou_jing uses the paper's signed polynomial kernel.",
                 call. = FALSE)
        return("zhou_jing")
    }
    if (method %in% c("chen_hall", "chen_hall_bartlett")) {
        if (is.null(kernel))
            return("biweight")
        return(match.arg(kernel, c(
            "biweight", "epanechnikov", "gaussian"
        )))
    }
    if (!is.null(kernel))
        stop("kernel is not used by this method.", call. = FALSE)
    NULL
}
