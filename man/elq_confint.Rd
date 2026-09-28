\name{elq_confint}
\alias{elq_confint}
\title{Quantile and Quantile-Difference Confidence Sets}
\description{Compute a closed confidence-set hull for one quantile, a quantile of
paired differences, or a difference of independent quantiles.}
\usage{
elq_confint(x, y = NULL, tau = .5, paired = FALSE, method = "el",
  conf.level = .95, interval = c("bonferroni", "profile"),
  bandwidth = NULL, kernel = NULL, adjustment = NULL,
  bartlett = "plugin", control = elq_control())
}
\arguments{
\item{x, y}{Numeric samples.}
\item{tau}{Quantile probability in (0,1).}
\item{paired}{Use the sample x-y. Equal lengths are required.}
\item{method}{One of the six methods returned by elq_methods().}
\item{conf.level}{Desired confidence level.}
\item{interval}{For independent differences: bonferroni subtracts marginal
confidence hulls at level 1-(1-conf.level)/2; profile directly inverts the
stepwise difference statistic and is available for el and ael.}
\item{bandwidth, kernel, adjustment, bartlett}{Method parameters; see elq_test.}
\item{control}{Settings from elq_control().}
}
\details{
Stepwise one-sample sets are calculated by count enumeration. Adimari inversion
uses the interpolated CDF. Zhou-Jing inversion enumerates the roots of the
piecewise cubic CDF analogue and permits disconnected sets. Chen-Hall inversion
uses its monotone branches around the zero of the smooth estimating equation.

The result is a closed interval hull; the one-sample components are retained
as an attribute. Closing endpoints or filling gaps can only enlarge the set.
There is no finite-sample coverage guarantee from chi-square calibration.
AEL can produce unbounded intervals. An empty set is returned as NA bounds
with an empty components table.

Bonferroni difference intervals are deliberately conservative and do not need to
agree with rejection by the profiled difference test.
}
\value{An elq_interval numeric vector with lower and upper bounds, confidence
level, method description, and, for one sample, the closed component intervals.}
\examples{
set.seed(9)
x <- rexp(40)
elq_confint(x, method = "zhou_jing")
elq_confint(x, rnorm(50), interval = "profile", method = "el")
}
