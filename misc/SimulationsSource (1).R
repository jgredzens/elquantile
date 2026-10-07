# ==============================================================================
#  Empirical Likelihood Inference for Quantiles — One-Sample / Paired Case
#
#  Methods:  NS-EL | Adjusted EL | Adimari | Zhou & Jing
#            Chen & Hall SEL | Bartlett BC-SEL
#            + Sign test | Wilcoxon | Bootstrap | Chen & Hall + Adjusted EL | Bartlett + Adjusted EL
#
#  References:
#    Owen (1988, 1990)               — Empirical likelihood foundations
#    Chen & Hall (1993)              — Smoothed EL for quantiles, Bartlett correction
#    Adimari (1998)                  — EL statistic for quantiles
#    Zhou & Jing (2003)              — Adjusted EL method for quantiles
#    Chen, Variyath & Abraham (2008) — Adjusted EL (AEL)
# ==============================================================================

library(ggplot2)
library(patchwork)

# ==============================================================================
#  SHARED CONSTANTS
# ==============================================================================

EL_METHODS <- c("standard", "adjusted", "adimari", "zhou_jing", "chen_hall", "bartlett")#, "chen_hall_adj", "bartlett_adj")

METHOD_LABELS <- c(
    # Baselines
    bootstrap   = "Bootstrap",
    # Classical non-parametric
    sign        = "Sign Test",
    wilcoxon    = "Wilcoxon",
    # Empirical likelihood
    standard    = "NS-EL",
    adjusted    = "Adjusted EL",
    adimari     = "Adimari",
    zhou_jing   = "Zhou & Jing",
    chen_hall   = "Chen & Hall",
    bartlett    = "Bartlett",
    chen_hall_adj = "Chen & Hall Adj",
    bartlett_adj = "Bartlett Adj"
)

# Colorblind-friendly palette (Wong 2011) + extras
METHOD_COLORS <- c(
    "Bootstrap"   = "#444444",  # Dark charcoal
    "Sign Test"   = "#2CA02C",  # Orange
    "Wilcoxon"    = "#56B4E9",  # Sky blue
    "NS-EL"       = "#D62728",  # Bluish green
    "Adjusted EL" = "#9467BD",  # Reddish purple
    "Adimari"     = "#8C564B",  # Vermillion
    "Zhou & Jing" = "#E377C2",  # Black
    "Chen & Hall" = "#7F7F7F",  # Deep blue
    "Bartlett"    = "#BCBD22",   # Yellow,
    "Chen & Hall Adj" = "#17BECF",
    "Bartlett Adj" = "#F2D0A4"
)

METHOD_ORDER <- c(
    "Bootstrap", "Sign Test", "Wilcoxon",
    "NS-EL", "Adjusted EL", "Adimari", "Zhou & Jing", "Chen & Hall", "Bartlett", "Chen & Hall Adj", "Bartlett Adj"
)

# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

# Silently return NA on error instead of stopping
safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

# Silverman's rule-of-thumb bandwidth (scale-robust)
silverman <- function(x) 0.9 * min(sd(x), IQR(x) / 1.34) * length(x)^(-1/5)

# NULL-safe fallback: use b if a is NA
`%||%` <- function(a, b) if (!is.na(a)) a else b

# Convert EL log-likelihood ratio to p-value
llr_to_pval <- function(llr) {
    if (is.na(llr) || !is.finite(llr) || llr < 0) return(NA_real_)
    pchisq(llr, df = 1, lower.tail = FALSE)
}

run_simulation <- function(sim_params = SIM, bw = BW) {
    # ... (rest of your existing run_simulation code) ...
}

# Helper function implementing the algorithm from the image
solve_lambda_newton <- function(g, max_iter = 100, eps = 1e-8) {
    n_star <- length(g)
    lambda <- 0  # Step 2: Initial value
    k <- 0
    gamma <- 1
    
    for (k in 1:max_iter) {
        # Step 3: Compute derivatives
        # R(lambda) = sum(log(1 + lambda * g))
        # First deriv (R_dot): sum( g / (1 + lambda * g) )
        # Second deriv (R_ddot): -sum( g^2 / (1 + lambda * g)^2 )
        
        arg <- 1 + lambda * g
        if (any(arg <= 0)) return(NA_real_) # Safety check for log domain
        
        R_dot <- sum(g / arg)
        R_ddot <- -sum(g^2 / arg^2)
        
        # Step 3 continued: Compute step Delta = -R_dot / R_ddot
        delta <- -R_dot / R_ddot
        
        # Convergence check
        if (abs(delta) < eps) return(lambda)
        
        # Step 4: Line search (Step size reduction)
        gamma <- 1
        # Check conditions: (1 + (lambda + gamma*delta)*g > 0) AND (R(new) > R(old))
        # Note: In Step 4 of the image, the goal is to maximize the dual, 
        # so R(new) must be > R(old).
        
        repeat {
            new_lambda <- lambda + gamma * delta
            new_arg <- 1 + new_lambda * g
            
            if (all(new_arg > 0)) {
                # Compare R values
                r_new <- sum(log(new_arg))
                r_old <- sum(log(arg))
                if (r_new > r_old) break
            }
            
            gamma <- gamma / 2
            if (gamma < 1e-10) return(NA_real_) # Line search failed to find improvement
        }
        
        # Step 5: Update
        lambda <- new_lambda
        # The image suggests gamma = (k+1)^-1/2, but Step 4's repeat loop
        # usually handles the specific step adjustment for convergence.
    }
    return(NA_real_) # Max iterations reached
}

solve_lambda_robust <- function(g) {
    lo <- 1 / (-max(g)) + 1e-8
    hi <- 1 / (-min(g)) - 1e-8
    obj <- function(l) sum(g / (1 + l * g))
    
    # 1. Try Newton-Raphson first (fast, usually best)
    lam <- tryCatch(solve_lambda_newton(g), error = function(e) NA_real_)
    if (!is.na(lam) && lam > lo && lam < hi) return(lam)
    
    # 2. Fall back to uniroot (guaranteed if bracketed)
    lam <- tryCatch(
        uniroot(obj, interval = c(lo, hi), tol = 1e-10)$root,
        error = function(e) NA_real_
    )
    if (!is.na(lam)) return(lam)
    
    # 3. Last resort: Brent over a tighter interval
    lam <- tryCatch(
        optimize(function(l) obj(l)^2, interval = c(lo, hi))$minimum,
        error = function(e) NA_real_
    )
    lam
}


# ==============================================================================
#  EL LOG-LIKELIHOOD RATIO
# ==============================================================================

#' Compute the empirical log-likelihood ratio for a quantile.
#'
#' @param x       Numeric vector of observations (differences for paired case).
#' @param theta   Hypothesised quantile value (theta_0).
#' @param t       Quantile level in (0, 1). Default 0.5 (median).
#' @param method  One of "standard", "adjusted", "adimari", "zhou_jing",
#'                "chen_hall", "bartlett", "chen_hall_adj", "bartlett_adj".
#' @param h       Bandwidth for smoothed methods. Defaults to silverman(x).
#' @return  Scalar -2 log R(theta), or NA if the optimization fails.
el_quantiles <- function(x, theta, t = 0.5, method = "standard", h = NULL, scale_ael = FALSE) {
    n <- length(x)
    if (is.null(h)) h <- silverman(x)

    # ── Closed-form methods (no lambda solve required) ─────────────────────────

    if (method == "adimari") {
        # Piecewise-linear interpolation of F*_n — Adimari (1998)
        x_s   <- sort(x)
        y_s   <- (2 * seq_len(n) - 1) / (2 * n)
        fn_star <- if (theta < x_s[1] | theta >= x_s[n])
            return(Inf)
        else
            approx(x_s, y_s, xout = theta)$y
        return(2 * n * (fn_star * log(fn_star / t) +
                        (1 - fn_star) * log((1 - fn_star) / (1 - t))))
    }

    if (method == "zhou_jing") {
        # Kernel-smoothed CDF — Zhou & Jing (2003)
        fn_hat <- mean(pnorm((theta - x) / h))
        if (fn_hat <= 0 || fn_hat >= 1) return(Inf)
        return(2 * n * (fn_hat * log(fn_hat / t) +
                        (1 - fn_hat) * log((1 - fn_hat) / (1 - t))))
    }

    # ── Methods solved via Lagrange multiplier (lambda) ────────────────────────
    # Step 1: estimating function
    # g <- switch(method,
    #             standard      = ,
    #             adjusted      = as.numeric(x <= theta) - t,
    #             chen_hall     = pnorm((theta-x) / h) - t,
    #             bartlett      = ,
    #             chen_hall_adj = pnorm((theta-x) / h) - t,
    #             bartlett_adj  = pnorm((theta-x) / h) - t,
    #             stop("Unknown method: ", method)
    # )
    # 
    # 
    # # AEL: append pseudo-observation — Chen, Variyath & Abraham (2008)
    # if (method %in% c("adjusted", "chen_hall_adj", "bartlett_adj")) {
    #   an <- max(1, log(n) / 2)
    #   g  <- c(g, -an * mean(g))
    # }
    # 
    # # Convex hull check: solution exists iff g straddles zero
    # if (max(g) <= 0 || min(g) >= 0) return(NA_real_)
    # 
    # # Solve dual equation: sum g_i / (1 + lambda * g_i) = 0
    # lambda <- safe(uniroot(
    #     f        = function(l) sum(g / (1 + l * g)),
    #     interval = c(1 / (-max(g)) + 1e-8, 1 / (-min(g)) - 1e-8),
    #     tol      = 1e-10
    # )$root)
    # 
    # if (is.na(lambda)) return(NA_real_)
    # 
    # d <- 1 + lambda * g
    # if (any(d <= 0)) return(NA_real_)
    # 
    # llr <- 2 * sum(log(d))
    # # # Rescale AEL back to n-observation scale 
    # if (scale_ael) {
    #     llr <- llr / (1 + an/n)
    # }
    
    
    # ── Methods solved via Lagrange multiplier (lambda) ────────────────────────
    g <- switch(method,
                standard      = as.numeric(x <= theta) - t,
                adjusted      = as.numeric(x <= theta) - t,
                chen_hall     = pnorm((theta-x) / h) - t,
                bartlett      = pnorm((theta-x) / h) - t,
                chen_hall_adj = pnorm((theta-x) / h) - t,
                bartlett_adj  = pnorm((theta-x) / h) - t,
                stop("Unknown method: ", method))
    
    # AEL: append pseudo-observation
    is_ael <- method %in% c("adjusted", "chen_hall_adj", "bartlett_adj")
    if (is_ael) {
        an <- max(1, log(n) / 2)
        g  <- c(g, -an * mean(g))
    }
    
    
    # --- UPDATED LOGIC ---
    # Use Newton-Raphson with line search for Adjusted methods
    # if (is_ael) {
    #     lambda <- solve_lambda_newton(g)
    # } else {
    #     # Standard methods can still use uniroot or Newton
    #     lambda <- tryCatch(
    #         uniroot(f = function(l) sum(g / (1 + l * g)),
    #                 interval = c(1/(-max(g)) + 1e-8, 1/(-min(g)) - 1e-8),
    #                 tol = 1e-10)$root,
    #         error = function(e) NA_real_
    #     )
    # }
    # 
    # 
    if (max(g) <= 0 || min(g) >= 0) return(NA_real_)
    lambda <- solve_lambda_robust(g)
    if (is.na(lambda)) return(NA_real_)
    
    
    # # 1. Convex hull failure
    # if (max(g) <= 0 || min(g) >= 0) {
    #     message("NA: convex hull — max(g)=", round(max(g),6), " min(g)=", round(min(g),6))
    #     return(NA_real_)
    # }
    # 
    # # 2. After solver
    # lambda <- solve_lambda_robust(g)
    # if (is.na(lambda)) {
    #     message("NA: solver failed — max(g)=", round(max(g),6), " min(g)=", round(min(g),6),
    #             " interval=(", round(1/(-max(g)),6), ",", round(1/(-min(g)),6), ")")
    #     return(NA_real_)
    # }
    # 
    # # 3. After d check
    # d <- 1 + lambda * g
    # if (any(d <= 0)) {
    #     message("NA: d<=0 after solve — lambda=", round(lambda,6))
    #     return(NA_real_)
    # }
    # ---------------------
    
    
    d <- 1 + lambda * g
    llr <- 2 * sum(log(d))
    
    if (scale_ael && is_ael) {
        an <- max(1, log(n) / 2)
        llr <- llr / (1 + an/n)
    }
    
    
    # === Bartlett correction — Chen & Hall (1993) ===
    # Correction evaluated at theta
    if (method %in% c("bartlett", "bartlett_adj")) {
        mu2      <- mean(g^2)
        if (mu2 < 1e-10) return(NA_real_)
        mu3      <- mean(g^3)
        mu4      <- mean(g^4)
        beta_hat <- (1 / 6) * (3 * mu4 / mu2^2 - 2 * mu3^2 / mu2^3)
        llr      <- llr / (1 + beta_hat / n)
    }
    llr
}

el_mean <- function(x, mu, ael=FALSE) {
    n      <- length(x)
    g      <- x - mu
    
    # AEL: append pseudo-observation
    if (ael) {
        a    <- max(1, log(n) / 2)   # Chen & Variyath's suggested a
        g    <- c(g, -a * mean(g))   # pseudo-observation shifts centroid
        n    <- n + 1L
    }
    
    if (max(g) <= 0 || min(g) >= 0) return(NA_real_)
    
    lambda <- solve_lambda_newton(g)
    
    if (is.na(lambda)) return(NA_real_)
    d <- 1 + lambda * g
    if (any(d <= 0)) return(NA_real_)
    2 * sum(log(d))
}
# el_mean(c(1,2,3), mu=2)
# ==============================================================================
#  CI BY GRID INVERSION
# ==============================================================================

#' Confidence interval for a quantile via EL ratio inversion.
#'
#' Evaluates el_quantiles on a coarse grid, then refines both boundary
#' crossings with uniroot.
#'
#' @param x,t,method,alpha,h  As in el_quantiles.
#' @param ngrid  Number of coarse grid points. Default 500.
#' @return  Numeric vector c(lower, upper), or c(NA, NA) on failure.
el_ci <- function(x, t = 0.5, method = "standard", alpha = 0.05,
                  h = NULL, ngrid = 500) {
    if (is.null(h)) h <- silverman(x)
    crit    <- qchisq(1 - alpha, df = 1)
    rng     <- diff(range(x))
    grid    <- seq(min(x) - 0.45 * rng, max(x) + 0.45 * rng, length.out = ngrid)
    eval_fn <- if (is.null(t))
        function(th) safe(el_mean(x, th))
    else
        function(th) safe(el_quantiles(x, th, t, method, h))
    vals <- sapply(grid, eval_fn)
    ok   <- which(!is.na(vals) & vals <= crit)
    if (!length(ok)) return(c(NA_real_, NA_real_))
    refine <- function(i1, i2) {
        safe(uniroot(
            f        = function(th) { v <- eval_fn(th); if (is.na(v)) NA_real_ else v - crit },
            interval = grid[c(i1, i2)],
            tol      = 1e-8
        )$root) %||% grid[round((i1 + i2) / 2)]
    }
    lo <- if (min(ok) > 1)     refine(max(1, min(ok) - 1), min(ok))     else grid[1]
    hi <- if (max(ok) < ngrid) refine(max(ok), min(ngrid, max(ok) + 1)) else grid[ngrid]
    c(lo, hi)
}

# ==============================================================================
#  SIGN TEST HELPERS
# ==============================================================================

#' Sign test p-value for quantile t.
#' Ties are excluded from the effective sample size (exact conditional test).
pval_sign <- function(x, theta0, t) {
    nn <- sum(x != theta0)
    if (nn == 0L) return(1)
    safe(binom.test(sum(x > theta0), nn, p = 1 - t)$p.value)
}

#' Sign test CI: order-statistic interval.
ci_sign <- function(x, t, alpha) {
    ds <- sort(x)
    n  <- length(x)
    ds[c(max(1L, qbinom(alpha / 2, n, t)),
         min(n,  qbinom(1 - alpha / 2, n, t)))]
}

# ==============================================================================
#  COMPARE ALL METHODS — NUMERIC TABLE
# ==============================================================================

#' Compute CIs and p-values for all methods on a single sample.
#'
#' @param x      Numeric vector of observations.
#' @param theta0 Null hypothesis value. Defaults to sample quantile.
#' @param t      Quantile level. Default 0.5.
#' @param alpha  Significance level. Default 0.05.
#' @param h      Bandwidth. Defaults to silverman(x).
#' @param ngrid  Grid size for EL CI inversion.
#' @param B      Bootstrap replicates.
#' @param seed   Optional RNG seed for reproducibility.
#' @return  data.frame with columns Method, P_Value, CI_Lower, CI_Upper,
#'          CI_Width, Reject.
compare_el_methods <- function(x, theta0 = NULL, t = 0.5, alpha = 0.05,
                                h = NULL, ngrid = 1000, B = 999, seed = NULL,
                               methods = EL_METHODS) {
    n    <- length(x)
    qhat <- as.numeric(quantile(x, t))
    if (is.null(theta0)) theta0 <- qhat
    if (is.null(h))      h      <- silverman(x)
    if (!is.null(seed))  set.seed(seed)

    el_q_methods <- intersect(methods, EL_METHODS)
    # ── EL ────────────────────────────────────────────────────────────────────
    llrs     <- sapply(el_q_methods, function(m) safe(el_quantiles(x, theta0, t, m, h)))
    pvals_el <- sapply(llrs, llr_to_pval)
    cis_el   <- setNames(
        lapply(el_q_methods, el_ci, x = x, t = t, alpha = alpha, h = h, ngrid = ngrid),
        el_q_methods
    )

    # ── Bootstrap ─────────────────────────────────────────────────────────────
    boot_q  <- replicate(B, quantile(sample(x, n, replace = TRUE), t))
    ci_boot <- as.numeric(quantile(boot_q, c(alpha / 2, 1 - alpha / 2)))
    pv_boot <- mean(abs(boot_q - mean(boot_q)) >= abs(qhat - theta0))

    # ── Classical non-parametric ───────────────────────────────────────────────
    wt <- if (t == 0.5)
        safe(wilcox.test(x, mu = theta0, conf.int = TRUE, conf.level = 1 - alpha))
    else NULL

    # ── Assemble ──────────────────────────────────────────────────────────────
    ci_list <- c(cis_el, list(
        bootstrap = ci_boot,
        sign      = ci_sign(x, t, alpha),
        wilcoxon  = if (!is.null(wt)) as.numeric(wt$conf.int) else c(NA_real_, NA_real_)
    ))
    pv_list <- c(as.list(pvals_el), list(
        bootstrap = pv_boot,
        sign      = pval_sign(x, theta0, t),
        wilcoxon  = if (!is.null(wt)) wt$p.value else NA_real_
    ))

    res <- data.frame(
        Method   = METHOD_LABELS[names(ci_list)],
        P_Value  = round(unlist(pv_list), 4),
        CI_Lower = round(sapply(ci_list, `[`, 1), 4),
        CI_Upper = round(sapply(ci_list, `[`, 2), 4),
        CI_Width = round(sapply(ci_list, function(ci) ci[2] - ci[1]), 4),
        Reject   = unlist(pv_list) < alpha,
        row.names = NULL
    )
    res$Method <- factor(res$Method, levels = rev(METHOD_ORDER))
    if (t != 0.5) res <- res[res$Method != "Wilcoxon", ]
    res
}

# ==============================================================================
#  PROFILE LIKELIHOOD PLOT
# ==============================================================================

#' Plot -2 log R(theta) vs theta for all EL methods.
#'
#' Overlays the chi-squared critical value and sign test CI boundaries.
#'
#' @param x      Numeric vector of observations.
#' @param t      Quantile level. Default 0.5.
#' @param alpha  Significance level. Default 0.05.
#' @param h      Bandwidth. Defaults to silverman(x).
plot_el_profile <- function(x, t = 0.5, alpha = 0.05, h = NULL) {
    if (is.null(h)) h <- silverman(x)
    crit <- qchisq(1 - alpha, df = 1)
    rng  <- diff(range(x))
    grid <- seq(min(x) - 0.3 * rng, max(x) + 0.3 * rng, length.out = 200)

    plot_data <- do.call(rbind, lapply(EL_METHODS, function(m) {
        data.frame(
            theta  = grid,
            llr    = sapply(grid, function(th) safe(el_quantiles(x, th, t, m, h))),
            method = METHOD_LABELS[m]
        )
    }))
    plot_data$method <- factor(plot_data$method, levels = METHOD_LABELS[EL_METHODS])

    # Sign test CI boundaries (order-statistic quantiles)
    n       <- length(x)
    ds      <- sort(x)
    sign_lo <- ds[max(1L, qbinom(alpha / 2,       n, t))]
    sign_hi <- ds[min(n,  qbinom(1 - alpha / 2,   n, t))]
    sign_df <- data.frame(
        xintercept = c(sign_lo, sign_hi),
        label      = c("Sign\nLower", "Sign\nUpper")
    )

    plot_colors <- setNames(
        METHOD_COLORS[METHOD_LABELS[EL_METHODS]],
        METHOD_LABELS[EL_METHODS]
    )

    ggplot(plot_data, aes(theta, llr, colour = method)) +
        geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
        geom_hline(yintercept = crit, linetype = "dashed",
                   colour = "red", linewidth = 0.8) +
        geom_vline(data = sign_df, aes(xintercept = xintercept),
                   linetype = "dotted", colour = "black", linewidth = 0.8) +
        annotate("text",
                 x = sign_df$xintercept, y = 9.5,
                 label = sign_df$label, size = 3,
                 hjust = c(1.1, -0.1), colour = "black") +
        scale_colour_manual(values = plot_colors, name = NULL) +
        labs(
            title    = "Profile empirical log-likelihood ratio",
            subtitle = sprintf(
                "Dashed red: \u03c7\u00b2(%d%%) cutoff  |  Dotted black: sign test %d%% CI",
                round((1 - alpha) * 100), round((1 - alpha) * 100)
            ),
            x = expression(theta),
            y = expression(-2 ~ log ~ R(theta))
        ) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
}

# ==============================================================================
#  DISTRIBUTION + CI COMPARISON PLOT
# ==============================================================================

#' KDE of the data with CIs for all methods, plus CI-width and p-value panels.
#'
#' @param x      Numeric vector of observations.
#' @param theta0 Null value. Default 0.
#' @param t      Quantile level. Default 0.5.
#' @param alpha  Significance level. Default 0.05.
#' @param h      Bandwidth. Defaults to silverman(x).
#' @param ngrid  Grid size for EL CI inversion.
#' @param B      Bootstrap replicates.
#' @param seed   Optional RNG seed.
#' @param label  Plot title prefix. Defaults to variable name.
#' @return  Returns plot list and the results data frame from compare_el_methods().
plot_el_comparison <- function(x, theta0 = NULL, t = 0.5, alpha = 0.05,
                               h = NULL, ngrid = 500, B = 999,
                               seed = NULL, label = NULL,
                               methods = c(EL_METHODS, "el_mean", "ttest"), nudge_p1_x = 0.75) {
    if (is.null(theta0)) theta0 <- 0
    if (is.null(h))      h      <- silverman(x)
    if (is.null(label))  label  <- deparse(substitute(x))

    # message("Computing CIs — this may take a few seconds...")
    n    <- length(x)
    qhat <- as.numeric(quantile(x, t))
    crit <- qchisq(1 - alpha, df = 1)

    # res  <- compare_el_methods(x, theta0 = theta0, t = t, alpha = alpha,
    #                            h = h, ngrid = ngrid, B = B, seed = seed)
    # res$Method <- factor(res$Method, levels = unique(res$Method))
    res_list <- list()
    
    el_q_methods <- intersect(methods, EL_METHODS)
    if (length(el_q_methods)) {
        res_q <- compare_el_methods(x, theta0 = theta0, t = t, alpha = alpha,
                                    h = h, ngrid = ngrid, B = B, seed = seed,
                                    methods = el_q_methods)
        res_q$Method <- factor(res_q$Method, levels = unique(res_q$Method))
        res_list[["el_quantiles"]] <- res_q
    }
    
    if ("el_mean" %in% methods) {
        ci_mean   <- el_ci(x, t = NULL, alpha = alpha, ngrid = ngrid)
        llr_mean  <- safe(el_mean(x, theta0))
        pval_mean <- if (!is.na(llr_mean)) 1 - pchisq(llr_mean, df = 1) else NA_real_
        res_list[["el_mean"]] <- data.frame(
            Method   = "EL Mean",
            CI_Lower = ci_mean[1],
            CI_Upper = ci_mean[2],
            CI_Width = ci_mean[2] - ci_mean[1],
            P_Value  = pval_mean,
            Reject   = !is.na(pval_mean) & pval_mean < alpha
        )
    }
    
    if ("ttest" %in% methods) {
        tt <- t.test(x, mu = theta0, conf.level = 1 - alpha)
        res_list[["ttest"]] <- data.frame(
            Method   = "t-test",
            CI_Lower = tt$conf.int[1],
            CI_Upper = tt$conf.int[2],
            CI_Width = diff(tt$conf.int),
            P_Value  = tt$p.value,
            Reject   = tt$p.value < alpha
        )
    }
    
    res <- do.call(rbind, res_list)
    rownames(res) <- NULL
    res$Method <- factor(res$Method, levels = unique(res$Method))
    
    # ── KDE ───────────────────────────────────────────────────────────────────
    bw_kde <- silverman(x)
    x_seq  <- seq(min(x) - 2 * sd(x), max(x) + 2 * sd(x), length.out = 300)
    kde_df <- data.frame(
        x = x_seq,
        y = sapply(x_seq, function(xi) mean(dnorm((xi - x) / bw_kde) / bw_kde))
    )
    y_max <- max(kde_df$y)

    # Label showing theta0 and qhat in left-to-right order
    if (theta0 <= qhat) {
        label_txt <- sprintf(
            "theta[0] == %.2f ~~ ',' ~~ hat(theta)[p] == %.2f",
            theta0, qhat
        )
    } else {
        label_txt <- sprintf(
            "hat(theta)[p] == %.2f ~~ ',' ~~ theta[0] == %.2f",
            qhat, theta0
        )
    }
    extra_colors <- c("EL Mean" = "#1F77B4", "t-test" = "#FF7F0E")
    all_colors   <- c(METHOD_COLORS, extra_colors)
    plot_colors  <- all_colors[levels(res$Method)]
    res$combo_label <- paste0(
        sprintf("%.3f", res$CI_Width),
        " | p ",
        ifelse(is.na(res$P_Value), "NA",
               ifelse(res$P_Value < 0.001, "< 0.001",
                      sprintf("= %.3f", res$P_Value))),
        ifelse(!is.na(res$Reject) & res$Reject == TRUE, " *", "")
    )
    res$p_label <- paste0(
        res$Method,
        " | p ",
        ifelse(is.na(res$P_Value), "NA",
               ifelse(res$P_Value < 0.001, "< 0.001",
                      sprintf("= %.3f", res$P_Value))),
        ifelse(!is.na(res$Reject) & res$Reject == TRUE, " *", "")
    )
    ci_df       <- res[!is.na(res$CI_Lower), ]
    ci_df$y_pos <- -rev(seq_len(nrow(ci_df))) * y_max * 0.13
    
    # ── Plot 1: KDE + CIs ─────────────────────────────────────────────────────
    p1 <- ggplot() +
        geom_ribbon(data = kde_df, aes(x = x, ymin = 0, ymax = y),
                    fill = "grey70", alpha = 0.2) +
        geom_line(data = kde_df, aes(x = x, y = y),
                  colour = "grey50", linewidth = 0.8) +
        geom_rug(data = data.frame(x = x), aes(x = x),
                 colour = "grey40", alpha = 0.5, linewidth = 0.3) +
        geom_segment(data = ci_df,
                     aes(x = CI_Lower, xend = CI_Upper,
                         y = y_pos, yend = y_pos, colour = Method),
                     linewidth = 1, lineend = "round") +
        geom_point(data = ci_df,
                   aes(x = CI_Lower, y = y_pos, colour = Method), size = 2) +
        geom_point(data = ci_df,
                   aes(x = CI_Upper, y = y_pos, colour = Method), size = 2) +
        geom_text(data = ci_df,
                  aes(x = CI_Upper, y = y_pos, label = p_label, colour = Method),
                  hjust = 0, nudge_x = nudge_p1_x, size = 5) +
        geom_vline(xintercept = theta0, linetype = "dashed",
                   colour = "#A32D2D", linewidth = 0.7) +
        geom_vline(xintercept = qhat, linetype = "dashed",
                   colour = "#185FA5", linewidth = 0.7) +
        # geom_label(
        #     aes(x = (theta0 + qhat) / 2, y = y_max * 1.05),
        #     label     = label_txt,
        #     parse     = TRUE,
        #     fill      = scales::alpha("white", 0.6),
        #     colour    = "black",
        #     size      = 3,
        #     label.size = 0.3,
        #     hjust     = 0.5
        # ) +
        scale_colour_manual(values = plot_colors) +
        scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
        labs(
            x        = expression(Delta),
            y        = "Density",
            title    = paste0(label, " KDE with ",
                              round((1 - alpha) * 100), "% CIs"),
            subtitle = bquote(
                p == .(t) ~ "|" ~ theta[0] == .(theta0) ~
                "|" ~ hat(theta)[0] == .(round(qhat,3)) ~ "|" ~ h == .(round(h, 3)) ~ "|" ~ n == .(n)
            )
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position  = "none",
            panel.grid.minor = element_blank(),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40")
        )

    # ── Plot 2: CI widths ─────────────────────────────────────────────────────
    p2 <- ggplot(res[!is.na(res$CI_Width), ],
                 aes(CI_Width, Method, fill = Method)) +
        geom_col(width = 0.6) +
        geom_text(aes(label = sprintf("%.3f", CI_Width)),
                  hjust = -0.1, size = 3, colour = "grey30") +
        scale_fill_manual(values = plot_colors) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
        labs(
            x     = expression("CI width  " ~ (theta[U] - theta[L])),
            y     = NULL,
            title = "CI width by method"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position    = "none",
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title         = element_text(size = 12, face = "bold")
        )

    # ── Plot 3: p-values ──────────────────────────────────────────────────────
    pval_df <- res[!is.na(res$P_Value), ]

    p3 <- ggplot(pval_df, aes(P_Value, Method, colour = Method)) +
        geom_vline(xintercept = alpha, linetype = "dashed",
                   colour = "grey40", linewidth = 0.6) +
        geom_segment(aes(x = 0, xend = P_Value, yend = Method),
                     linewidth = 0.7, alpha = 0.4) +
        geom_point(aes(shape = Reject), size = 3) +
        geom_text(
            aes(label = ifelse(
                is.na(P_Value), "",
                ifelse(P_Value < 0.001, "<0.001", sprintf("%.3f", P_Value))
            )),
            hjust = -0.3, size = 3, colour = "grey30"
        ) +
        scale_colour_manual(values = plot_colors) +
        scale_shape_manual(
            values = c("TRUE" = 17, "FALSE" = 16),
            labels = c("TRUE" = "Reject", "FALSE" = "Fail to reject")
        ) +
        scale_x_continuous(
            limits = c(0, max(pval_df$P_Value, na.rm = TRUE) * 1.3),
            expand = expansion(mult = c(0, 0.05))
        ) +
        annotate("text",
                 x = alpha, y = 0,
                 label  = sprintf("alpha == %.2f", alpha),
                 parse  = TRUE,
                 colour = "grey40",
                 size   = 3,
                 hjust  = -0.1,
                 vjust  = 0) +
        labs(
            x     = expression(italic(p) * "-value"),
            y     = NULL,
            title = bquote("p-values at" ~ theta[0] == .(theta0)),
            shape = NULL
        ) +
        expand_limits(y = -0.5) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position    = "none",
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title         = element_text(size = 12, face = "bold")
        )

    # ── Plot 4: CI width + p-value combined ───────────────────────────────────
    
    p4 <- ggplot(res, aes(CI_Width, Method, fill = Method)) +
        geom_col(width = 0.6) +
        geom_text(
            aes(label = combo_label, colour = Reject),
            hjust = -0.1, size = 5
        ) +
        scale_fill_manual(values = plot_colors) +
        scale_colour_manual(
            values = c("TRUE" = "firebrick", "FALSE" = "grey30"),
            guide  = "none"
        ) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
        labs(
            # x        = expression("CI width  " ~ (theta[U] - theta[L])),
            x = NULL,
            y        = NULL,
            title    = "CI width by method",
            subtitle = bquote("(*) significant at" ~ alpha == .(alpha))
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position    = "none",
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title         = element_text(size = 12, face = "bold"),
            plot.subtitle      = element_text(size = 10, colour = "grey40")
        )
    
    
    # ── Combine ───────────────────────────────────────────────────────────────
    combined <- 
        p1 / (p2 | p3) +
            plot_layout(heights = c(2, 1.2)) +
            plot_annotation(
                caption = bquote(
                    h == .(round(h, 3)) ~ "|" ~
                    B == .(B) ~ "|" ~
                    "grid" == .(ngrid) ~ "|" ~
                    chi[.(paste0(round((1 - alpha) * 100), "%"))]^2 * "(1)" == .(round(crit, 3))
                ),
                theme = theme(plot.caption = element_text(size = 9, colour = "grey50"))
            )

    # invisible(res)
    list(
        combined = combined,
        p1 = p1, 
        p2 = p2,
        p3 = p3,
        p4 = p4,
        res_df = res
    )
    
}
