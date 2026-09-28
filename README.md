# elquantile

Empirical likelihood quantile inference, power, and sample-size planning in R.

This is a modular research package whose inference functions depend only on
standard R. Building its optional vignette uses knitr, rmarkdown, and Pandoc.
It implements the six quantile EL constructions discussed in the supplied
papers and supports one sample, paired data, two independent samples, and
several independent groups. The method engines, profiling, intervals,
post-hoc inference, and planning calculations live in separate modules so
the package can be extended.

## Install

From the project directory:

~~~sh
R CMD INSTALL elquantile
~~~

Alternatively, install the enclosed source archive:

~~~r
install.packages("elquantile_0.1.0.tar.gz", repos = NULL, type = "source")
library(elquantile)
~~~

The archive is in the project's validation directory. R >= 3.6 is declared;
the enclosed checks were run on R 4.3.3 on Linux. No CRAN submission has been
made. Replace the placeholder author/maintainer metadata before release.

## Tutorial vignette

The executable tutorial covers the inference workflows, both Bartlett
corrections, post-hoc inference, power, sample size, and numerical diagnostics.
After installing the built source archive, open it with:

~~~r
vignette("quantile-inference", package = "elquantile")
~~~

Its source is vignettes/quantile-inference.Rmd, with references in
vignettes/references.bib. To build it from the package's parent directory,
first install knitr and rmarkdown in R and make Pandoc available:

~~~r
install.packages(c("knitr", "rmarkdown"))
~~~

~~~sh
R CMD build elquantile
R CMD INSTALL elquantile_0.1.0.tar.gz
~~~

Installing the raw source directory alone does not build the vignette.
The project also includes a rendered HTML copy in documentation/.

## Methods

| Method argument | Construction | Default smoothing |
|---|---|---|
| el | Ordinary binary estimating-equation EL | None |
| adimari | Midpoint-interpolated CDF in the explicit Bernoulli KL statistic | No bandwidth |
| chen_hall | Smoothed estimating equations, solving the scalar EL dual | Biweight, scale × n^(-1/2) |
| chen_hall_bartlett | Chen-Hall divided by 1 + Bartlett factor/n | Biweight, scale × n^(-3/4) |
| zhou_jing | Signed polynomial kernel CDF substituted into Bernoulli KL | Original signed kernel, scale × n^(-1/2) |
| ael | Add estimating pseudo-observation -a × mean(g), a=max(1,log(n)/2) | None |

The default scale is IQR/(2*qnorm(0.75)). Pass a bandwidth explicitly to
reproduce a paper's unscaled rule. Chen-Hall also supports ordinary
Epanechnikov on [-1,1] and Gaussian kernels. The original Chen-Hall numerical
example uses a rescaled Epanechnikov kernel on [-sqrt(5),sqrt(5)]; reproduce
that by multiplying its reported bandwidth by sqrt(5) with our Epanechnikov.

The two meanings of “adjusted” are separate methods: Zhou-Jing is not
pseudo-observation AEL. See inst/doc/METHODS.md for equations and source mapping.

## One sample

~~~r
library(elquantile)
set.seed(27)
x <- rexp(80)

# Is the population median log(2)?
elq_one_sample(x, q0 = log(2), tau = 0.5,
               method = "chen_hall_bartlett", conf.int = TRUE)

# A confidence set's closed hull; individual components are retained.
ci <- elq_confint(x, tau = 0.9, method = "zhou_jing")
ci
attr(ci, "components")

# Exact binomial benchmark for a population quantile.
elq_sign_test(x, q0 = log(2), tau = 0.5)
~~~

The Bartlett option defaults to the empirical-moment correction. Use
bartlett="limit" for the simpler probability-only factor. They are different
corrections with different theoretical accuracy statements.

## Two samples: choose the scientific target

~~~r
before <- rnorm(70)
after <- before + rexp(70) - log(2)

# Quantile of within-person differences: q_tau(before - after).
elq_two_sample(before, after, paired = TRUE, delta = 0,
               method = "adimari", conf.int = TRUE)

# Difference of two independent population quantiles: q_tau(A) - q_tau(B).
A <- rnorm(60)
B <- rnorm(90, mean = 0.4, sd = 2)
elq_two_sample(A, B, delta = 0, tau = 0.5, method = "ael")
elq_confint(A, B, method = "ael", interval = "profile")
~~~

The quantile of a difference is generally not the difference of quantiles.
For paired inference, the rows are paired before any sorting.

All methods support independent-difference confidence intervals through a
Bonferroni difference of marginal confidence hulls. Direct profile-difference
intervals are currently implemented for el and ael. These alternatives are
labelled in the returned interval; Bonferroni intervals need not invert the
profile test decision.

## Several groups and post-hoc inference

~~~r
groups <- list(A = rnorm(60), B = rnorm(80, sd = 2),
               C = rnorm(70, mean = 0.5, sd = 3))
elq_k_sample(groups, tau = 0.5, method = "chen_hall")

# One Holm family across ALL pairs and supplied quantile levels.
elq_posthoc(groups, tau = c(0.25, 0.5, 0.75),
             method = "el", p.adjust.method = "holm", conf.int = TRUE)

# Formula interface.
d <- data.frame(value = unlist(groups),
                treatment = rep(names(groups), lengths(groups)))
elq_anova(value ~ treatment, data = d, tau = 0.5, method = "ael")
~~~

These compare quantiles, not means. Unequal group distributions and unequal
sample sizes are allowed under the stated regularity conditions. Groups must
be independent. This version does not fit crossed factorial effects,
repeated-measures models or quantile regression.

Post-hoc intervals use separate Bonferroni constructions. They are not
inversions of Holm-adjusted p-values. Post-hoc testing is not gated by an
omnibus rejection.

## Analytical power and sample size

Power requires a population quantile and the density at that quantile for
each group. A raw quantile difference alone is insufficient.

~~~r
# Normal groups with medians 0 and 0.5, standard deviations 1 and 2.
q <- c(0, 0.5)
f <- c(dnorm(0), dnorm(0)/2)
elq_power(n = c(100, 150), quantiles = q, density = f)

# Minimum integer base size for a 1:2 allocation and 80% approximate power.
plan <- elq_sample_size(power = 0.8, quantiles = q, density = f,
                        ratio = c(1, 2))
plan
plan$previous.power

# Three groups, testing equality of their medians.
elq_power(n = c(80, 100, 120), quantiles = c(0, 0.2, 0.6),
           density = c(dnorm(0), dnorm(0)/2, dnorm(0)))

# Single-sample upper-tail test at a nonmedian quantile.
elq_sample_size(power = 0.9,
  quantiles = qnorm(0.9) + 0.3, density = dnorm(qnorm(0.9)),
  tau = 0.9, null = qnorm(0.9), alternative = "greater")
~~~

For paired planning, provide the quantile and density of the **paired
differences** as one sample. The returned n is the number of pairs.
For example, normally distributed differences with median 0.4 and standard
deviation 1 have density dnorm(0) at their median.

The analytical calculator uses local asymptotic normal/noncentral chi-square
theory. It is shared across the regular EL variants and does not claim
method-specific small-sample power. An optional contrast matrix supports
more general planning questions; arbitrary-contrast EL testing is not yet
implemented.

elq_density(x, tau) supplies a pilot KDE estimate when needed. Its uncertainty
is not automatically propagated into planning. Assess plausible density
values and verify designs by simulation.

## Exact count-based power

Ordinary EL, AEL, and the exact sign test depend only on the number of
observations below a one-sample null quantile. Their rejection probability
can therefore be enumerated without simulation.

~~~r
# Alternative X ~ N(0.5, 1), testing median q0=0:
p.alt <- pnorm(0, mean = 0.5)
elq_power_exact(n = 40, p = p.alt, method = "el")
elq_power_exact(n = 40, p = p.alt, method = "ael")
elq_sample_size_exact(power = 0.8, p = p.alt,
                      method = "sign", max.n = 500)
~~~

Here p means F_alternative(q0), not tau or the alternative quantile.
The result reports both power and actual null size. Exact enumeration does
not make asymptotic EL cutoffs exact: a rare-quantile EL test may exceed
nominal size. Set max.null.size in exact sample-size planning if required.
The search scans integers because discrete power can oscillate.

## Method-specific simulation planning

~~~r
gens <- list(function(n) rnorm(n),
              function(n) rnorm(n, mean = 0.5, sd = 2))
sim <- elq_power_sim(c(80, 120), generators = gens,
                     method = "zhou_jing", nsim = 1000, seed = 431)
sim
sim$mc.interval
sim$failure.bounds

grid <- elq_sample_size_sim(c(40, 80, 120), generators = gens,
  ratio = c(1, 2), power = 0.8, method = "el",
  nsim = 1000, seed = 432)
grid$table

paired.generator <- function(n) {
  baseline <- rnorm(n)
  cbind(baseline + rnorm(n, 0.4, 1), baseline)
}
elq_power_sim(60, paired.generator, paired = TRUE,
               method = "chen_hall", nsim = 500, seed = 433)
~~~

A generator satisfying the null estimates size; an alternative generator
estimates power. In asymmetric distributions, center by the target population
quantile rather than the mean. Use the same seed and generators to compare
the deterministic inference methods on the same datasets.

Simulation output includes individual outcomes, warnings, failed replicates,
Monte Carlo standard error, a Wilson interval, and failure bounds. Power is
reported conditional on successful replicates; failures are never silently
discarded. With a supplied seed, the caller's RNG state is restored.

Selection from a simulation grid is exploratory. It is not a guarantee of
minimal sample size or achieved population power. Confirm the selected size
using an independent larger run. Pointwise Monte Carlo confidence limits
are not simultaneous guarantees across a searched grid.

## Likelihood curves and lower-level access

~~~r
curve <- elq_profile(x, theta = seq(0.1, 2, length.out = 200),
                     tau = 0.5, method = "adimari")
plot(curve)
elq_methods()
elq_dual(c(-0.5, 0.5), weights = c(12, 18))
~~~

## Numerical and statistical scope

- Continuous data, independent sampling, unique quantiles and positive local
  densities are required. Ties, missing and nonfinite observations are
  rejected; no observations are dropped automatically.
- Ordinary EL returns infinity for an infeasible constraint. This is a valid
  statistic outcome, not an optimizer failure.
- AEL adjusts estimating equations; it does not invent a measured observation.
  Its confidence set can be unbounded.
- The Zhou-Jing kernel is signed. Its CDF analogue can be nonmonotone and
  confidence sets can be disconnected. Values outside the open probability
  branch receive an infinite statistic; no clipping hides this issue.
- Stepwise group profiles are enumerated exactly. Smooth group profiles use
  knots and numerical refinement. Inspect refinement.gap; stability does not
  certify a global optimum.
- Group extensions have first-order calibration. The original one-sample
  higher-order guarantees do not automatically transfer to profiling,
  multiplicity, random bandwidths or arbitrary custom settings.
- A “Bartlett” method label is not a universal guarantee of second-order
  accuracy. Kernel order, bandwidth rate and smoothness matter.

## Package structure and extension

| Module | Responsibility |
|---|---|
| R/control.R | Validation, controls, RNG handling |
| R/methods.R | Method metadata, kernels and bandwidth rules |
| R/dual.R | Scalar EL solver, binary EL/AEL and Bartlett factors |
| R/model.R | Prepared one-group criterion objects |
| R/profiling.R | Independent-group constrained minimization |
| R/intervals.R | One-sample sets and difference intervals |
| R/inference.R | User APIs and post-hoc inference |
| R/power.R | Analytical and exact planning |
| R/simulation.R | Method-specific simulation planning |

See inst/doc/EXTENDING.md for the engine contract and a change checklist.
The source project includes validation scripts and checks. To run:

~~~sh
R CMD build elquantile
R CMD check --no-manual elquantile_0.1.0.tar.gz
~~~

The project studies directory contains executable validation simulations.
These are implementation diagnostics, not evidence that one method is
uniformly superior or that the package establishes a novel statistical method.
