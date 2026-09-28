\name{elq_test}
\alias{elq_test}
\alias{elq_one_sample}
\alias{elq_two_sample}
\alias{elq_k_sample}
\alias{elq_anova}
\alias{elq_posthoc}
\alias{elq_sign_test}
\title{Empirical Likelihood Tests for Population Quantiles}
\description{One-sample, paired-difference, independent two-sample and independent
k-sample quantile inference using six empirical likelihood constructions.
An exact binomial sign test is also supplied as a separate benchmark.}
\usage{
elq_test(x, y = NULL, group = NULL, tau = .5, null = 0,
  paired = FALSE, method = "el",
  alternative = c("two.sided", "less", "greater"),
  conf.int = FALSE, conf.level = .95, bandwidth = NULL, kernel = NULL,
  adjustment = NULL, bartlett = "plugin", control = elq_control())
elq_one_sample(x, q0 = 0, ...)
elq_two_sample(x, y, delta = 0, paired = FALSE, ...)
elq_k_sample(x, group = NULL, ...)
elq_anova(formula, data, ...)
elq_posthoc(x, group = NULL, tau = .5, method = "el",
  p.adjust.method = c("holm", "bonferroni", "none"),
  conf.int = FALSE, conf.level = .95, bandwidth = NULL,
  kernel = NULL, adjustment = NULL, bartlett = "plugin",
  control = elq_control())
elq_sign_test(x, q0 = 0, tau = .5,
  alternative = c("two.sided", "less", "greater"), conf.level = .95)
}
\arguments{
\item{x, y}{Numeric samples; x may instead be a named list of independent samples.
Ties within an analyzed sample, missing values and nonfinite values are rejected.}
\item{group}{Grouping labels when x is a numeric response vector.}
\item{tau}{Prespecified probability in (0,1). A vector of distinct probabilities
is accepted by elq_posthoc; all comparisons form one multiplicity family.}
\item{null, q0, delta}{Null quantile for one sample or paired differences; null
quantile difference for two independent groups. Use zero for k-group equality.}
\item{paired}{Analyze the quantile of x-y, preserving observation pairing.
This is not the difference of marginal quantiles.}
\item{method}{One of el, adimari, chen_hall, chen_hall_bartlett, zhou_jing, ael.}
\item{alternative}{Two-sided, lower or upper alternative. Directional alternatives
are available only for a scalar target and use signed-root normal calibration.}
\item{conf.int}{Return confidence intervals. Intervals are always two-sided.
Use post-hoc intervals rather than an omnibus interval with more than two groups.}
\item{conf.level}{Confidence level. In post-hoc inference this is the simultaneous level.}
\item{bandwidth}{NULL for the documented method-specific rule, or a positive
scalar or vector in group order. Only smoothed methods use bandwidths.}
\item{kernel}{NULL for the method default. Chen-Hall allows biweight,
epanechnikov or gaussian. Zhou-Jing uses its signed polynomial kernel.}
\item{adjustment}{Positive pseudo-observation factor for ael, scalar or per group.
NULL uses max(1,log(n)/2).}
\item{bartlett}{For chen_hall_bartlett: plugin for empirical moments at the usual
sample quantile, or limit for the probability-only limiting factor.}
\item{control}{Settings from elq_control().}
\item{formula}{A formula response ~ group with one grouping variable.}
\item{data}{Data frame containing the formula variables.}
\item{p.adjust.method}{Holm, Bonferroni or no multiplicity adjustment.}
\item{...}{Other arguments passed to elq_test or elq_k_sample.}
}
\details{
The one-sample statistic is the method's pointwise negative twice log relative
likelihood, or its explicit KL analogue. For independent groups the statistic
is the common-quantile constrained sum minus the unrestricted sum of minima.
For two groups, x is shifted by the hypothesized quantile difference.
The reference degrees of freedom are one for a scalar target and k-1 for
k-group equality.

The finite-sample unrestricted minimum of a stepwise EL objective need not be
zero. The group comparison retains this baseline. The empirical CDF cells
are enumerated exactly for el and ael. Adimari uses convex minimization on
every interpolation segment. Chen-Hall and Zhou-Jing use numerical profiling;
returned refinement gaps are diagnostics, not global optimality certificates.

Ordinary EL assigns an infinite statistic to an infeasible moment constraint.
AEL adds one estimating pseudo-observation, not a new measured data value.
The signed Zhou-Jing CDF analogue can leave [0,1]; outside its open real-valued
KL branch the statistic is set to infinity and the CDF is never clipped.

All groupwise extensions use first-order calibration under independent,
continuous populations with unique quantiles, positive local densities, a
fixed group count and nonvanishing sample fractions. Smoothing needs further
kernel, density and bandwidth assumptions. The one-sample higher-order
theorems in the papers do not establish higher-order accuracy for these group
extensions. These are quantile comparisons, not classical mean ANOVA.

Post-hoc tests profile each pair separately, including under partial nulls.
Holm adjustment is across all supplied pairs and quantiles. Optional intervals
are separate Bonferroni differences of marginal confidence hulls, not inversions
of Holm decisions. An omnibus rejection is not required before post-hoc tests.

The sign benchmark uses exact binomial probabilities. Its order-statistic
interval is a two-sided equal-tail interval and may be unbounded. An observation
equal to q0 is rejected as an ambiguous boundary case.
}
\value{Tests return an elq_test object inheriting from htest, with statistic,
p.value, estimate, target, calibration, group sizes, bandwidths and numerical
diagnostics. Post-hoc inference returns a data frame with raw and adjusted
p-values, and optional confidence bounds.}
\references{
Adimari (1998), \doi{10.1080/00949659808811873}.
Chen and Hall (1993), \doi{10.1214/aos/1176349256}.
Zhou and Jing (2003), Annals of the Institute of Statistical Mathematics 55, 689--703.
Chen, Variyath and Abraham (2008), Adjusted empirical likelihood and its properties.
}
\seealso{\code{\link{elq_methods}}, \code{\link{elq_confint}}, \code{\link{elq_power}}}
\examples{
set.seed(18)
x <- rnorm(40)
y <- rnorm(45, sd = 2)
elq_one_sample(x, q0 = 0, method = "adimari", conf.int = TRUE)
elq_two_sample(x, y, method = "el")
elq_k_sample(list(A = x, B = y, C = rnorm(50, sd = 3)))
elq_posthoc(list(A = x, B = y), conf.int = TRUE)
}
