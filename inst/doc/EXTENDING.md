# Extending elquantile

The public APIs depend on a prepared one-group model, not directly on a
particular solver. Keep statistical definitions separate from optimization
and from simulation.

## Add a method

1. Add its metadata row in R/methods.R. Use an unambiguous identifier and
   cite the precise source equation or label it an experimental extension.
2. Add kernel or bandwidth helpers there if needed. Kernels must document
   their support, order, normalization and whether they are signed.
3. Add a branch to .eq_model in R/model.R. A model stores the sample, n,
   tau, method, sample.quantile, and a vectorized value(theta) function.
   It must also supply:
   - minimum: the unrestricted criterion infimum;
   - center: a minimizer or suitable central root;
   - bounds and knots for profiling;
   - costs for a stepwise method, indexed by count+1.
4. Route profiling in R/profiling.R. Do not apply smooth optimization to
   a step function. Do not call a finite grid globally exact. Expose numerical
   diagnostics and infeasibility explicitly.
5. Add confidence-set inversion in R/intervals.R. Preserve disconnected
   components when relevant, or explicitly return a conservative hull.
6. Add mathematical identity and regression tests, including affine
   transformation, boundary and infeasibility cases. Check paired equivalence
   and two-group/list interfaces through the public API.
7. Update the README, Rd pages, method mapping and release notes; rerun
   R CMD build and R CMD check.

The registry is deliberately internal and explicit. Adding a row alone is
not enough to register a statistically valid method. An optional public
custom-method registration interface can be added later once its invariants
and calibration contract are settled.

## Extend study planning

The analytical engine takes quantiles, densities and a contrast matrix.
It can support future contrast-based inference without changing the planning
formula. A new method with a different limit distribution must not silently
reuse the noncentral chi-square approximation.

The simulation engine uses user-supplied generators and the public test API.
It can assess new methods once the inference route is implemented. Preserve
failure counts, uncertainty, seed behavior and reproducible raw outputs.

Potential next additions include joint multi-quantile EL, general linear
contrasts, bootstrap calibration, clustered data, covariate adjustment and
bootstrap simultaneous intervals. These require new statistical definitions
and validation; they are not implemented by merely adding an option.

## Compatibility and release

Keep the existing meanings of tau, null, paired and the sign of a difference.
For paired data, n always counts pairs. For independent designs n contains
one sample size per group. Confidence intervals and p-value adjustments
must retain their explicit construction labels.

The package is version 0.1.0 research software. Author/maintainer fields are
placeholders, not attribution to the attached papers' authors. Supply real
metadata and a release repository before distributing a maintained package.
No user data or attached paper PDFs are bundled.
