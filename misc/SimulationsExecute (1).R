# ==============================================================================
#  Simulation driver — coverage probability across methods and distributions
#
#  Generates paired differences, computes p-values for all EL and classical
#  methods, and summarises empirical coverage at multiple nominal levels.
#
#  Usage:
#    1. Adjust SIM and BW constants at the bottom of this file.
#    2. Run the full script; results are saved to sim_results.rds.
#    3. Load saved results with: readRDS("sim_results.rds")
# ==============================================================================
# getwd()
# 
# 
library(ggplot2)
library(dplyr)

sign_test_p <- function(x, mu = 0) {
    x_clean <- x[x != mu]          # remove exact ties
    binom.test(sum(x_clean > mu), length(x_clean))$p.value
}

theme_springer <- function(base_size = 17) {
    theme_bw(base_size = base_size) %+replace%
        theme(
            text = element_text(size=base_size),
            # --- Panel ---
            panel.background  = element_rect(fill = "white", colour = NA),
            plot.background   = element_rect(fill = "white", colour = NA),
            panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.4),
            panel.grid.major  = element_line(colour = "grey85", linewidth = 0.2),
            panel.grid.minor  = element_blank(),
            
            # --- Axes ---
            axis.ticks        = element_line(colour = "black", linewidth = 0.3),
            axis.ticks.length = unit(0.12, "cm"),
            axis.text         = element_text(size = rel(0.9), colour = "black"),
            axis.title        = element_text(size = rel(1.0), colour = "black"),
            axis.title.x      = element_text(margin = margin(t = 6)),
            axis.title.y      = element_text(margin = margin(r = 6), angle = 90),
            
            # # --- Legend ---
            # legend.background = element_rect(fill = "white", colour = "grey70",
            #                                  linewidth = 0.3),
            # legend.key        = element_rect(fill = "white", colour = NA),
            # legend.key.size   = unit(0.8, "lines"),
            # legend.text       = element_text(size = rel(0.9)),
            # legend.title      = element_text(size = rel(0.9)),
            # legend.margin     = margin(3, 5, 3, 5),
            legend.position   = "right",
            # legend.position   = "none",
            
            # --- Facets ---
            strip.background  = element_rect(fill = "grey92", colour = "black",
                                             linewidth = 0.3),
            strip.text        = element_text(size = rel(0.9)),
            
            # --- Titles ---
            plot.title        = element_text(size = rel(1.05), hjust = 0.5),
            plot.subtitle     = element_text(size = rel(0.9), hjust = 0.5,
                                             colour = "grey30"),
            plot.caption      = element_text(size = rel(0.8), hjust = 0,
                                             colour = "grey40"),
            plot.margin       = margin(6, 6, 6, 6),
            
            complete = TRUE
        )
}

{
  source("../Downloads/SimulationsSource (1).R")
  set.seed(13)
  x <- rexp(50)
  # x <- rnorm(50)
  x <- x - median(x)
  # x <- x - mean(x)
  # x <- c(x - mean(x), rep(100,3))
  plot_el_profile(x, t = 0.5) + ylim(c(-.5,8)) + xlim(c(-.5,1)) + theme_springer() + labs(title=NULL, subtitle = NULL) + theme(
      legend.text = element_text(size=20)
  )
  ggsave("ProfileEL.eps", width = 420, height = 160, units = "mm", device = cairo_ps)
  results <- plot_el_comparison(x, theta0 = 0, t = 0.5,
                     label = "Centered exp(1)", h=length(x)^(-0.5),
                     # methods = c("ttest", "el_mean", "NS-EL", "bartlett")
                     )

  print(
      results$p1 + labs(x=NULL, title=NULL)  + theme_springer() | 
          results$p4+ labs(x=NULL, title=NULL) + theme_springer() + theme(panel.grid.major = element_blank(),
                                                                          panel.grid.minor = element_blank(),
                                                                          panel.background = element_blank())
  )
  ggsave("example_exp.eps", width = 420, height = 160, units = "mm", device = cairo_ps)
  # ggsave("example_norm.eps", width = 420, height = 160, units = "mm", device = cairo_ps)
  # t.test(x)
  # # Sample data
  # emplik::el.test(x, mu=0)#$Pval
}
{
    # Country GDP example: https://wid.world/data/
    # https://cran.r-project.org/web/packages/wid/wid.pdf
    wealth <- read.csv("../Downloads/national-wealth.csv", sep=";")
    w_data <- wealth %>% mutate(
        d = v2024 - v2014, 
    ) %>%
        filter(!grepl("Russian Federation|Belarus|Ukraine|Italy|United Kingdom|Moldova", Country)) %>%
        arrange(d)
    # theta_0 <- median(w_data$d)
    theta_0 <- 40000
    w_data %>% dplyr::summarise(
        n = n(),
        shapiro_p = shapiro.test(d)$p.value,
        mean = mean(d),
        median = median(d),
        t_p = t.test(d, mu=theta_0)$p.value,
        elmean_p = emplik::el.test(d, mu=theta_0)$Pval,
        wilcox_p = wilcox.test(d, mu = theta_0)$p.value,
        sign_p = sign_test_p(d, mu=theta_0),
        el_q25_p = llr_to_pval(el_quantiles(d, theta = theta_0, t=0.25, method="bartlett", h=n^(-1/2))),
        el_q50_p = llr_to_pval(el_quantiles(d, theta = theta_0, t=0.5, method="bartlett", h=n^(-1/2))),
        el_q75_p = llr_to_pval(el_quantiles(d, theta = theta_0, t=0.75, method="bartlett", h=n^(-1/2)))
    ) %>%
    arrange(-t_p) %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
    # ggplot(w_data) +
    #     geom_point(aes(x = v2014, y=v2024)) +
    #     geom_abline()
    # plot_el_comparison(x=w_data$delta,
    #                    theta0 = median(w_data$delta),
    #                    # theta0 = 0,
    #                    t = 0.5,
    #                    h=0.5*length(w_data$delta)^(-0.2),
    #                    label = "2024 vs 2014 individual wealth data by country \u0394")
    # plot_el_profile(w_data$delta, t = 0.5)+ ylim(c(-.1,50)) + xlim(c(0,100000))
    plot_methods <- c("ttest", "el_mean", "standard", "adjusted", "bartlett")
    results <- plot_el_comparison(
        x=w_data$d,
        # theta0 = median(w_data$delta),
        theta0 = theta_0,
        t = 0.5,
        h=length(w_data$d)^(-1/2),
        methods = plot_methods,
        label = "2024 vs 2014 individual wealth data by country \u0394", nudge_p1_x = 10000)
    # print(results$p1 + results$p4)# + plot_layout(widths=c(2,1)))
    results$p1 + labs(title=NULL) + theme_springer() + theme(legend.position = "none") 
    # ggsave("country_deltas_withfilter_p1.eps", width = 200, height = 160, units = "mm", device = cairo_ps)
    # summary(wealth)
    # median(wealth$v2024) - median(wealth$v2014)
    # mean(w_data$delta)
    # t.test(w_data$delta, mu=80000)
    # emplik::el.test(w_data$delta, mu=80000)$Pval
    # myfun <- function(theta, x) {
    #     emplik::el.test(x, mu = theta)
    # }
    # 
    # # findUL finds both upper and lower CI bounds automatically
    # # result <- 
    # emplik::findUL(step = 20, fun = myfun, MLE = mean(w_data$delta-mean(w_data$delta)), x = w_data$delta - mean(w_data$delta))
    # result$Low   # lower bound
    # result$Up    # upper bound
    # t-test can't be used
    # EL means is still usable
}

# ==============================================================================
#  QUICK EXPLORATORY EXAMPLE  (Sleep dataset, n = 10)
# ==============================================================================

{
    data(sleep)
    delta_sleep <- sleep$extra[sleep$group == 2] - sleep$extra[sleep$group == 1]

    cat("\n=== Sleep Study (n =", length(delta_sleep), ") ===\n")
    cat("Sample median difference:", round(median(delta_sleep), 3), "\n")
    cat("Silverman h:             ", round(silverman(delta_sleep), 3), "\n\n")

    # Numeric comparison table
    print(compare_el_methods(delta_sleep, theta0 = 1.5, t = 0.5, seed = 1))

    # Profile likelihood
    plot_el_profile(delta_sleep, t = 0.5) + ylim(c(-.1,8)) + xlim(c(-0.5,5))

    # KDE + CI + p-value panels
    plot_el_comparison(delta_sleep, theta0 = 0, t = 0.5,
                       label = "Sleep \u0394", seed = 1)
}


# ==============================================================================
#  SIMULATION INFRASTRUCTURE
# ==============================================================================

# ── Method groups ─────────────────────────────────────────────────────────────

METHODS_NONSMOOTH <- c("standard", "adjusted", "adimari")
METHODS_SMOOTH    <- c("chen_hall", "bartlett", "zhou_jing", "chen_hall_adj", "bartlett_adj")

# ── Supported distributions ───────────────────────────────────────────────────
# Each entry provides:
#   rgen — random sample generator: function(n) -> numeric vector
#   qfun — quantile function:       function(t) -> scalar true quantile

DISTRIBUTIONS <- list(
    normal = list(
        rgen = function(n) rnorm(n),
        qfun = function(t) qnorm(t)
    ),
    chisq3 = list(
        rgen = function(n) rchisq(n, df = 3),
        qfun = function(t) qchisq(t, df = 3)
    ),
    exp1 = list(
        rgen = function(n) rexp(n),
        qfun = function(t) qexp(t)
    )
)

# ── Simulation parameters ─────────────────────────────────────────────────────
# For final paper results set nsim = 10000 and uncomment all n / t / alpha values.

SIM <- list(
    dist_names = names(DISTRIBUTIONS),          # distributions to iterate over
    n_vec      = c(20, 50, 100),
    t_vec      = c(0.05, 0.1, 0.5, 0.9, 0.95),
    alpha_vec  = c(0.05, 0.1),
    nsim       = 10000,    # increase to 10000 for paper results
    seed       = 42
)

# ── Bandwidth strategies ──────────────────────────────────────────────────────
# "n^{-1/5}" and "n^{-1/3}" are scale-free rates; include at least one
# scale-aware option (Silverman) for each run.

BW <- list(
    # "silverman" = function(x) silverman(x),
    # "0.5*n^{-0.2}"  = function(x) 0.5*length(x)^(-1/5),
    # "n^{-0.2}"  = function(x) length(x)^(-1/5),
    # "1.5n^{-0.2}"  = function(x) 1.5*length(x)^(-1/5),
    "2n^{-0.2}"  = function(x) 2*length(x)^(-1/5)#,
    # "0.5*n^{-0.75}"  = function(x) 0.5*length(x)^(-3/4),
    # "1n^{-0.75}"  = function(x) length(x)^(-3/4),
    # "1.5n^{-0.75}"  = function(x) 1.5*length(x)^(-3/4),
    # "2n^{-0.75}"  = function(x) 2*length(x)^(-3/4)#,
    # "n^{-0.5}"  = function(x) length(x)^(-0.5),
    # "n^{-0.75}"  = function(x) length(x)^(-0.75),
    # "n^{-1}"  = function(x) length(x)^(-1)#,
    # "0.3"       = function(x) 0.3
)

# ==============================================================================
#  ONE SIMULATION ITERATION
# ==============================================================================

#' Compute p-values for all methods on a single sample.
#'
#' @param x       Numeric vector (one simulated sample).
#' @param theta0  True quantile value under H0.
#' @param t       Quantile level.
#' @param h_vals  Named list of bandwidth functions (name -> function(x)).
#' @return  Named list of p-values, one per method/bandwidth combination.
one_iter_pval <- function(x, theta0, t, h_vals) {
    results <- list()

    # Non-smoothed EL methods (no bandwidth needed)
    for (m in METHODS_NONSMOOTH) {
        results[[m]] <- llr_to_pval(
            safe(el_quantiles(x, theta = theta0, t = t, method = m))
        )
    }
    # results[["scaled_ael"]] <- llr_to_pval(
    #     safe(el_quantiles(x, theta = theta0, t = t, method = "adjusted", scale_ael = TRUE))
    # )
    # Smoothed EL methods — one result per bandwidth strategy
    for (m in METHODS_SMOOTH) {
        for (bw_name in names(h_vals)) {
            h   <- h_vals[[bw_name]](x)
            key <- paste0(m, ".", bw_name)
            results[[key]] <- llr_to_pval(
                safe(el_quantiles(x, theta = theta0, t = t, method = m, h = h))
            )
        }
    }
    
    # Classical non-parametric
    results[["sign_test"]] <- pval_sign(x, theta0, t)
    if (t == 0.5) {
        results[["wilcoxon"]] <- safe(wilcox.test(x - theta0)$p.value)
    }

    results
}

# ==============================================================================
#  SIMULATION LOOP
# ==============================================================================

#' Run the full coverage simulation.
#'
#' Iterates over all combinations of distribution × n × t, draws nsim samples,
#' computes p-values, and summarises empirical coverage at each alpha level.
#'
#' @param sim_params  List with fields: dist_names, n_vec, t_vec, alpha_vec,
#'                    nsim, seed.
#' @param bw          Named list of bandwidth functions (passed to one_iter_pval).
#' @return  data.frame with columns: dist, n, t, gamma, key, coverage.
run_simulation <- function(sim_params = SIM, bw = BW) {
    set.seed(sim_params$seed)
    
    conditions <- expand.grid(
        dist = sim_params$dist_names,
        n    = sim_params$n_vec,
        t    = sim_params$t_vec,
        stringsAsFactors = FALSE
    )
    
    all_results <- vector("list", nrow(conditions))
    
    for (cond_idx in seq_len(nrow(conditions))) {
        dist_name <- conditions$dist[cond_idx]
        n_val     <- conditions$n[cond_idx]
        t_val     <- conditions$t[cond_idx]
        
        dist_obj  <- DISTRIBUTIONS[[dist_name]]
        theta0    <- dist_obj$qfun(t_val)   
        
        smooth_keys    <- as.vector(outer(METHODS_SMOOTH, names(bw), paste, sep = "."))
        nonsmooth_keys <- c(METHODS_NONSMOOTH)
        classical_keys <- if (t_val == 0.5) c("sign_test", "wilcoxon") else "sign_test"
        all_keys       <- c(nonsmooth_keys, smooth_keys, classical_keys)
        
        pval_mat <- matrix(
            NA_real_,
            nrow     = sim_params$nsim,
            ncol     = length(all_keys),
            dimnames = list(NULL, all_keys)
        )
        
        for (i in seq_len(sim_params$nsim)) {
            x    <- dist_obj$rgen(n_val)
            iter <- one_iter_pval(x, theta0, t_val, bw)
            for (key in all_keys) pval_mat[i, key] <- iter[[key]]
        }
        
        # --- MODIFICATION START ---
        # Calculate how many NAs occurred for each key once per condition block
        na_counts <- colSums(is.na(pval_mat))
        # --- MODIFICATION END ---
        
        cond_rows <- lapply(sim_params$alpha_vec, function(alpha) {
            data.frame(
                dist      = dist_name,
                n         = n_val,
                t         = t_val,
                gamma     = 1 - alpha,
                key       = all_keys,
                coverage  = 1 - colMeans(pval_mat < alpha, na.rm = TRUE),
                na_count  = na_counts, # Added column here
                row.names = NULL
            )
        })
        all_results[[cond_idx]] <- do.call(rbind, cond_rows)
        
        cat(sprintf(
            "[%d/%d] dist=%-8s n=%2d t=%.2f done\n",
            cond_idx, nrow(conditions), dist_name, n_val, t_val
        ))
    }
    
    do.call(rbind, all_results)
}

# ==============================================================================
#  RUN AND SAVE
# ==============================================================================

{
    sim_results <- run_simulation(SIM, BW)
    saveRDS(sim_results, "sim_results_v14.rds")
}
#     result <- tidyr::pivot_wider(dplyr::distinct(sim_results), names_from = c("t", "n"), values_from = "coverage", values_fn = mean)
write.csv(sim_results, "sim_results_14.csv", row.names = F)
# result <- readRDS("sim_results_v6.rds")
    

sim_results %>%
    select(-coverage) %>%
    # filter(!(t %in% c(0.05, 0.95))) %>%
    tidyr::pivot_wider(., names_from = c("t", "key"), values_from = "na_count") %>% 
    write.table("clipboard", sep="\t", dec=".", row.names = FALSE)
sim_results %>%
    select(-na_count) %>%
    # filter(!(t %in% c(0.05, 0.95))) %>%
    tidyr::pivot_wider(., names_from = c("t", "key"), values_from = "coverage") %>% 
    write.table("clipboard", sep="\t", dec=".", row.names = FALSE)
