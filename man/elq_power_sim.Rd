\name{elq_power_sim}
\alias{elq_power_sim}
\alias{elq_sample_size_sim}
\title{Simulation Power and Sample-Size Grid Assessment}
\description{Estimate method-specific rejection probabilities with retained
failures and Monte Carlo uncertainty.}
\usage{
elq_power_sim(n, generators, tau = .5, null = 0, paired = FALSE,
  method = "el", nsim = 1000L, alpha = .05,
  alternative = c("two.sided", "less", "greater"),
  seed = NULL, keep = TRUE, ...)
elq_sample_size_sim(n.grid, generators, power = .8, ratio = NULL,
  paired = FALSE, criterion = c("estimate", "lower"),
  nsim = 1000L, seed = NULL, ...)
}
\arguments{
\item{n}{Group sizes or the number of pairs.}
\item{generators}{For one sample, a function(n) returning observations.
For independent groups, a list of such functions. For paired data, one
function(n) returning an n by 2 numeric matrix.}
\item{tau}{Quantile probability.}
\item{null}{Null quantile or quantile difference. Use zero for k-group equality.}
\item{paired}{Analyze the quantile of generated paired differences.}
\item{method}{EL method identifier.}
\item{nsim}{Number of independent simulation replicates, at least two.}
\item{alpha}{Nominal significance level.}
\item{alternative}{Two-sided or directional scalar alternative.}
\item{seed}{Optional nonnegative integer. When supplied, the caller's RNG state
is restored on exit. NULL consumes the current RNG stream.}
\item{keep}{Retain individual outcomes, errors, warnings and refinement gaps.}
\item{n.grid}{Distinct integer base sizes to assess.}
\item{power}{Requested power for selection from the evaluated grid.}
\item{ratio}{Allocation ratios, default equal allocation.}
\item{criterion}{Select using the point estimate or the lower pointwise Wilson limit.}
\item{...}{Additional inference parameters for elq_power_sim, or simulation
parameters for elq_sample_size_sim.}
}
\details{
The result estimates Type I error when the generating distributions satisfy
the null, and power otherwise. There is no automatic recentering: construct
generators with the intended population quantiles.

Failed replicates are recorded. The reported estimate uses valid replicates;
failure bounds additionally report the result if all failures were nonrejections
or rejections. An infinite EL statistic caused by an infeasible constraint is a
valid outcome with asymptotic p-value zero, not a silently removed error.

Sample-size simulation evaluates the entire supplied grid. Selection requires
zero failed replicates at that candidate. Pointwise Monte Carlo intervals
do not provide simultaneous or post-selection guarantees. Confirm a selected
design in an independent larger run.
}
\value{elq_power_sim returns an elq_power object with power, MCSE, a Wilson
interval, failure count, worst-case bounds and optional replicate data.
elq_sample_size_sim returns the selected group sizes (NULL if none), a summary
table, and all simulation objects.}
\examples{
elq_power_sim(30, function(n) rnorm(n, .5), nsim = 20, seed = 12)
elq_power_sim(c(30, 40),
  list(function(n) rnorm(n), function(n) rnorm(n, .5)),
  nsim = 10, seed = 13)
}
