\name{elquantile-package}
\alias{elquantile-package}
\alias{elquantile}
\docType{package}
\title{Empirical Likelihood Quantile Inference and Study Planning}
\description{
Inference for one quantile, a quantile of paired differences, independent
quantile differences, and equality of independent-group quantiles. The package
also provides post-hoc comparisons and analytical, exact-count, and simulation
power and sample-size calculations.
}
\details{
Start with the executable tutorial:
\code{vignette("quantile-inference", package = "elquantile")}.
It explains scientific targets, method choices, Bartlett corrections,
confidence-set conventions, diagnostics, and planning assumptions.

Inference functions use standard R. The optional \pkg{knitr} and
\pkg{rmarkdown} packages and Pandoc are needed to build the vignette,
not to run the inference functions. Install a built source archive to
include the rendered vignette.

Continuous-population assumptions and method-specific regularity conditions
apply. The current input validator rejects ties and missing or nonfinite
values. Independent-group corrections use first-order calibration; one-sample
higher-order guarantees are not asserted for these extensions.
}
\seealso{
\code{\link{elq_methods}}, \code{\link{elq_test}},
\code{\link{elq_confint}}, \code{\link{elq_posthoc}},
\code{\link{elq_power}}, \code{\link{elq_sample_size}},
\code{\link{elq_power_exact}}, \code{\link{elq_power_sim}}.
}
\keyword{package}
