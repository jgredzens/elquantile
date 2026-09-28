\name{elq_statistic}
\alias{elq_statistic}
\alias{elq_profile}
\alias{elq_methods}
\alias{elq_bandwidth}
\alias{elq_kernel}
\alias{elq_bartlett}
\alias{elq_dual}
\alias{elq_control}
\alias{elq_density}
\alias{plot.elq_profile}
\title{Quantile Likelihood Engines, Kernels and Numerical Controls}
\description{Lower-level building blocks for the inference and study-planning APIs.}
\usage{
elq_statistic(x, theta, tau = .5, method = "el", bandwidth = NULL,
  kernel = NULL, adjustment = NULL, bartlett = "plugin",
  control = elq_control())
elq_profile(x, theta, tau = .5, method = "el", bandwidth = NULL,
  kernel = NULL, adjustment = NULL, bartlett = "plugin",
  control = elq_control())
elq_methods()
elq_bandwidth(x, method = "chen_hall", scale = NULL)
elq_kernel(u, kernel = c("biweight", "epanechnikov", "gaussian", "zhou_jing"),
  cumulative = FALSE)
elq_bartlett(x, tau = .5, bandwidth = NULL, kernel = "biweight",
  type = c("plugin", "limit"))
elq_dual(g, weights = rep(1, length(g)), adjustment = NULL, tolerance = 1e-10)
elq_control(grid.points = 81L, refine = TRUE, tolerance = 1e-9,
  max.cells = 2e6, warn = TRUE)
elq_density(x, tau = .5, bandwidth = NULL)
\method{plot}{elq_profile}(x, conf.level = .95, ...)
}
\arguments{
\item{x}{Numeric observations, or an elq_profile object for plot.}
\item{theta}{Finite vector of proposed quantile values.}
\item{tau}{Quantile probability.}
\item{method}{Method identifier; see elq_methods().}
\item{bandwidth}{Positive smoothing bandwidth or NULL.}
\item{kernel}{Kernel identifier. Chen-Hall defaults to biweight. Zhou-Jing
requires its signed kernel.}
\item{adjustment}{Positive factor a for the pseudo-observation -a times the
weighted mean estimating value. NULL uses no adjustment in elq_dual.}
\item{bartlett, type}{plugin or limit for the Bartlett factor.}
\item{control}{An elq_control object.}
\item{scale}{Positive scale multiplying the sample-size bandwidth rule;
NULL uses IQR divided by twice the standard-normal 0.75 quantile.}
\item{u}{Kernel arguments; infinite values are allowed.}
\item{cumulative}{Return the integrated kernel instead of its derivative.}
\item{g}{Finite scalar estimating values.}
\item{weights}{Nonnegative frequency weights, one per estimating value.}
\item{tolerance}{Numerical root/minimization tolerance.}
\item{grid.points}{Initial profile search grid size, at least 21.}
\item{refine}{Evaluate a finer profile grid and report its objective gap.}
\item{max.cells}{Limit on exact confidence-set cells or profile knots.}
\item{warn}{Emit tail-sparsity and numerical-refinement warnings.}
\item{conf.level}{Chi-square reference level displayed by the plot.}
\item{...}{Other arguments passed to the plotting function.}
}
\details{
The scalar dual normalizes estimating values before root finding and searches
within the positive-denominator region. Returned probabilities are aggregate
masses when frequency weights are used; an adjustment adds a final mass.
No parameter-space extension is performed for an infeasible ordinary EL problem.

For Zhou-Jing the kernel is (A*u^2+B) on [-1,1], with
A=(21-9*sqrt(21))/8 and B=(-3+3*sqrt(21))/8. It is signed; its integral is not
a monotone CDF. Replacing it by an ordinary Gaussian kernel does not reproduce
the paper's construction or assumptions.

Default bandwidths equal scale*n^(-1/2) for Chen-Hall and Zhou-Jing and
scale*n^(-3/4) for Bartlett-corrected Chen-Hall. These are transparent starting
rules, not estimated optimal bandwidths. The papers' higher-order results
have additional smoothness and bandwidth conditions. In particular, h=n^(-1)
does not satisfy n*h/log(n) tending to infinity.

The plug-in Bartlett factor is m4/(2*m2^2)-m3^2/(3*m2^3), with moments of
G((qhat-X)/h)-tau evaluated at the usual type-1 sample quantile. The limiting
factor is (1-tau+tau^2)/(6*tau*(1-tau)). The statistic is divided by 1+factor/n.
These two factors have different one-sample accuracy guarantees; no generic
second-order guarantee is made for k-group profiles.

elq_density uses a Gaussian KDE at the type-1 sample quantile. Its default
bandwidth uses the same IQR scale times n^(-1/5). It is a pilot planning input,
not an uncertainty-adjusted density estimate.
}
\value{elq_statistic returns statistics; elq_profile returns a data frame;
elq_methods returns method metadata; elq_bandwidth and elq_bartlett return
scalars; elq_kernel returns a numeric vector; elq_dual returns the statistic,
multiplier, masses, feasibility and moment residual; elq_control returns settings;
elq_density returns an estimate with diagnostic attributes.}
\examples{
set.seed(3)
x <- rnorm(50)
elq_methods()
elq_statistic(x, c(-.2, 0, .2), method = "ael")
elq_bartlett(x)
elq_dual(c(-1, 2), weights = c(5, 3))
}
