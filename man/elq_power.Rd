\name{elq_power}
\alias{elq_power}
\alias{elq_sample_size}
\alias{elq_power_exact}
\alias{elq_sample_size_exact}
\title{Analytical Power and Sample-Size Planning for Quantile Inference}
\description{Local asymptotic power for one or more independent quantiles,
and exact one-sample count enumeration for unsmoothed EL, AEL or the sign test.}
\usage{
elq_power(n, quantiles, density, tau = .5, contrast = NULL, null = 0,
  alpha = .05, alternative = c("two.sided", "less", "greater"))
elq_sample_size(power = .8, quantiles, density,
  ratio = rep(1, length(quantiles)), tau = .5, contrast = NULL,
  null = 0, alpha = .05,
  alternative = c("two.sided", "less", "greater"),
  min.n = 2L, max.n = 1000000L)
elq_power_exact(n, p, tau = .5, method = c("el", "ael", "sign"),
  alpha = .05, alternative = c("two.sided", "less", "greater"),
  adjustment = NULL)
elq_sample_size_exact(power = .8, p, tau = .5,
  method = c("el", "ael", "sign"), alpha = .05,
  alternative = c("two.sided", "less", "greater"),
  min.n = 2L, max.n = 2000L, adjustment = NULL,
  max.null.size = NULL)
}
\arguments{
\item{n}{Integer sample sizes at least two. For analytical power, a scalar
is replicated across groups; exact enumeration uses a single n.}
\item{quantiles}{Alternative population quantile in each independent group,
in measurement units. For paired inference supply the quantile of differences.}
\item{density}{Positive population density at each alternative quantile.
A scalar is replicated. For paired inference use the density of differences.}
\item{tau}{Common quantile probability.}
\item{contrast}{Full-row-rank contrast matrix. NULL tests a single quantile
against null, or all group quantiles against the last group. General contrasts
are planning calculations; the inference API implements equality and pair differences.}
\item{null}{Null contrast value, scalar or one per row of contrast.}
\item{alpha}{Nominal significance level.}
\item{alternative}{Two-sided or directional; directional power requires one contrast.}
\item{power}{Requested power for sample-size planning.}
\item{ratio}{Positive allocation ratios. Group sizes are ceiling(base*ratio/min(ratio)).}
\item{min.n, max.n}{Integer search bounds for the base size or one-sample n.}
\item{p}{Alternative probability P(X <= q0), not the quantile value or an effect size.}
\item{method}{Count-based test whose rejection region is enumerated.}
\item{adjustment}{Pseudo-observation factor for AEL; NULL uses max(1,log(n)/2).}
\item{max.null.size}{Optional upper bound on the enumerated null rejection
probability when selecting an exact-count sample size.}
}
\details{
With v=tau*(1-tau), the asymptotic covariance is diagonal with entries
v/(n_j*f_j(q_j)^2). If d=C*q-null, the noncentrality is
d'*(C*V*C')^(-1)*d. Two-sided power is the upper-tail probability of a
noncentral chi-square with rank(C) degrees of freedom at the central cutoff.
Directional power uses a signed normal approximation.

This local-alternative approximation is common to the regular EL variants.
It is not a finite-sample or Bartlett-specific power formula. It requires
densities; quantile differences alone do not determine power. Large effects,
extreme quantiles and density estimation uncertainty need simulation assessment.
Paired planning is one-sample planning for the distribution of differences,
so n is the number of pairs, not twice that number.

Analytical sample size uses monotonicity along the allocation path and verifies
the preceding integer base size. Exact count planning scans every n because
discrete power need not be monotone.

The exact functions enumerate M~Binomial(n,p). They report actual null size
using p=tau as well as power. This is exact enumeration of the specified test;
EL and AEL still use asymptotic cutoffs and can have null rejection probability
different from alpha. For a finite-sample level guarantee use the sign test.
Enumeration is limited to n<=100000.
}
\value{Power functions return an elq_power object. Planning functions return
elq_sample_size objects with group sizes, achieved power and search diagnostics.
Exact power additionally returns the full count rejection region and null size.}
\examples{
elq_power(c(100, 120), quantiles = c(0, .5), density = dnorm(0))
elq_sample_size(quantiles = c(0, .5), density = dnorm(0), ratio = c(1, 2))
elq_power_exact(40, p = pnorm(0, mean = .5), method = "el")
elq_sample_size_exact(power = .8, p = .3, method = "sign", max.n = 200)
}
