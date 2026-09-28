# Statistical definitions and source mapping

Write tau for the target probability and q for a candidate quantile.
All scalar engines return a nonnegative statistic on the scale -2 log R.
One-sample tests use the raw statistic; independent-group tests subtract the
unrestricted sum of group minima before chi-square calibration.

## Ordinary EL

g_i(q) = I(X_i <= q) - tau. For m observations at or below q and 0<m<n:

W(q) = 2 [m log{m/(n tau)} + (n-m) log{(n-m)/(n(1-tau))}].

At m=0 or n the moment constraint is infeasible, so W=Inf. The finite
algebraic endpoint limit of Bernoulli KL is not the empirical likelihood
value at an infeasible constraint. At noninteger n*tau, the unrestricted
minimum over q need not be zero.

## Adimari (1998)

Source: supplied article, printed pages 87-88, equations (2.2)-(2.3).
DOI: https://doi.org/10.1080/00949659808811873

Inside [X_(1),X_(n)), linearly interpolate the values (i-1/2)/n at the
ordered observations. Substitute this F* into 2n D(F* || tau).
The implementation retains the original feasible domain and gives an
infinite statistic outside it. Confidence sets are returned as closed hulls.

The article establishes the chi-square limit. Its statement that coverage
error should be O(n^-1) is explicitly a conjecture near the end of Section 3;
the package does not upgrade it to a proved general rate.

## Chen-Hall (1993) and Bartlett corrections

Primary source checked in addition to the attachments:
https://www.songxichen.com/Uploads/Files/Publication/Chen-Hall-93-AS.pdf
DOI: https://doi.org/10.1214/aos/1176349256

Set g_i(q)=G((q-X_i)/h)-tau, solve
sum g_i/(1+lambda*g_i)=0 with positive denominators, then return
2 sum log(1+lambda*g_i).

The default biweight density is 15(1-u^2)^2/16 on [-1,1]; its integrated
kernel is 1/2 + (15/16)(u - 2u^3/3 + u^5/5) on that interval.
This is the kernel used for the Chen-Hall comparator in Zhou-Jing's Section 3.

The full plug-in Bartlett factor, from Chen-Hall Section 4 (printed p.1172),
is b_hat = mu4_hat/(2 mu2_hat^2) - mu3_hat^2/(3 mu2_hat^3).
Moments are averages of [G((q_hat-X_i)/h)-tau]^j at the ordinary sample
quantile. The implementation uses R's type-1 empirical quantile.
Divide W by 1+b_hat/n. The simpler limiting correction is
b0=(1-tau+tau^2)/(6 tau(1-tau)); it is available with bartlett="limit".

The paper distinguishes their coverage rates. Its full higher-order result
requires the stated density smoothness, kernel order r, n*h/log(n) -> Inf,
and bounded n^3*h^(2r). For a second-order compact kernel, h proportional
to n^(-3/4) is a permitted rate. The uncorrected O(n^-1) theorem uses bounded
n*h^r, suggesting n^(-1/2) for r=2. Our data-scaled defaults are practical
starting rules; we do not assert every random bandwidth inherits the
deterministic-bandwidth higher-order theorem.

The ordinary Epanechnikov option is supported on [-1,1]. Chen-Hall's
Section 5 uses an equivalent rescaled kernel on [-sqrt(5),sqrt(5)].
Convert their h to sqrt(5)*h when using our Epanechnikov option.

## Zhou-Jing (2003)

Source: supplied article, equations (2.8)-(2.9), printed pp.691-692.
Title: Adjusted empirical likelihood method for quantiles.
Annals of the Institute of Statistical Mathematics 55, 689-703.

Use F_h(q)=mean G((q-X_i)/h) and W=2n D(F_h(q)||tau).
This substitutes the smooth CDF into the binary formula; it does not solve
the Chen-Hall estimating-equation dual.

The paper's kernel is K(u)=(A u^2+B) I(|u|<=1),
A=(21-9 sqrt(21))/8, B=(-3+3 sqrt(21))/8.
On [-1,1], G(u)=1/2+B u+A u^3/3. It is signed and its integral can
overshoot [0,1]. In particular, a conventional Gaussian CDF is not a faithful
substitute for this kernel's special cancellation condition.

We implement the open real-valued branch 0<F_h(q)<1; invalid probability
values, including exact endpoints, receive W=Inf. This boundary convention
is explicit because the displayed paper formula contains logarithms and
does not provide an extrapolation beyond [0,1]. We never silently truncate
F_h. This convention can matter for small samples and extreme quantiles.

The CDF analogue is piecewise cubic between X_i +/- h. Its level roots
are enumerated with polynomial roots for one-sample confidence sets, so
disconnected components are retained before returning their closed hull.

The paper's condition is n*h/log(n) -> Inf and bounded n*h^2, with additional
smoothness and kernel assumptions. It is not n*h*log(n) -> Inf.
The h=n^(-1) setting explored in the ICoMS manuscript does not satisfy
the first condition. It is available through an explicit bandwidth but
does not carry that higher-order theorem.

## Chen-Variyath-Abraham (2008) AEL

Source: supplied “Adjusted empirical likelihood and its properties”,
Section 3, equations (5)-(6), printed manuscript pp.6-8.

Append g_(n+1)=-a_n*mean(g_i), a_n=max(1,log(n)/2), then solve the same
EL dual with n+1 estimating values. The denominator normalization is n+1
and the pseudo-observation is included in the log-ratio sum.

This is different from Zhou-Jing's use of the word adjusted. For quantile
inference the original g_i are binary; frequency aggregation preserves the
exact objective. The pseudo-observation gives finite values even when all
original signs agree. It need not make a small-sample quantile procedure
exactly calibrated.

## The supplied ICoMS manuscript

Gredzens, Valeinis and Akopjana: Advantages of Empirical Likelihood Over
Wilcoxon Test for Paired Two-Sample Inference.

Its paired problem is implemented by forming X-Y and applying the one-sample
engine. Quantiles of differences and differences of marginal quantiles are
separate targets. A Wilcoxon signed-rank test is not added as an unrestricted
quantile engine: median interpretation needs symmetry. Use stats::wilcox.test
as an explicitly labelled benchmark when those assumptions are appropriate.

The package uses normalized relative likelihood prod(n*p_i), or
prod((n+1)*p_i) for AEL, and positive -2 log R statistics. It follows these
definitions rather than copying inconsistent normalization/sign expressions
from a draft.

## Independent groups: package extensions

For a proposed common quantile t, sum group criteria. Subtract the sum
of their separate unrestricted minima. For q_A-q_B=delta, compare
group A at t+delta to group B at t. Under standard first-order regularity
the reference distribution has k-1 degrees of freedom.

These independent-group profiles are package extensions of the supplied
one-sample constructions. The attachments do not themselves prove a
multi-group higher-order Bartlett theorem. Exact cell enumeration is used
for stepwise methods. Numerical smooth profiles return refinement diagnostics.
The method names identify engines, not a promise of identical higher-order
theory for every design.

## Power and sample size

Let V_j=tau(1-tau)/(n_j f_j(q_j)^2). For a full-row-rank contrast C and
d=C*q-null, the local noncentrality is d' [C diag(V_j) C']^(-1) d.
The noncentral chi-square upper tail gives approximate two-sided power.
For two groups this reduces to the squared quantile difference divided by
V_1+V_2. A one-sample or paired design uses a single variance.

First-order analytical power is shared across methods. Exact count-based
power is available for one-sample EL/AEL/sign because the statistic depends
only on M~Binomial(n,F_alt(q0)). Method-specific simulation can assess all
six methods and designs. No power calculation can determine a finite-sample
alternative distribution from a quantile difference alone.
