# IdSS (development version)

- Organize source files by method and consistently format R code and examples.
- Remove unused internal experiments and the unnecessary stringr dependency.
- Replace the expanded Gaussian-mixture derivatives by posterior-weighted
  component scores and curvature; preserve the public interface.
- Share inference selection in recursive SVARs and use explicit type predicates.
- Correct the sign-restriction help example layout and consolidate release notes.

# IdSS 0.3.1

Scientific review of the Structural Shocks companion course:
- Correct LP normalization to a fixed impact innovation, use usable-sample
  HAC scaling, and expose an optional HAC bandwidth for LP and LP-IV.
- Correct sign/zero samplers: centered spherical draws, Haar QR signs, all
  bootstrap replications, null-space ranks and separate inequality checks.
  Return each shock's own summaries and report empty admissible sets clearly.
- Retain early IRFs when maximizing FEV over windows beginning after step 1.
- Reject invalid Monte Carlo covariance draws and complete the requested draws.
- Use mean bootstrap bias, actual intercept vectors and linear stability
  shrinkage for the bootstrap-after-bootstrap procedure.
- Start plotted response horizons at impact zero; clarify single-draw FEVD
  dimensions and use relevant instruments in the TSLS help example.
- Clarify international macroeconomic data units.
- No datasets, exports or new dependencies are added.

# IdSS 0.3.0

Adds ten documented public functions for latent-state models:
`ms_stationary`, `ms_simulate`, `ms_filter` (Hamilton filter and backward
smoother), `fit_msar` (two-state constrained EM/ECM), `ms_forecast_moments`,
`ms_forecast_draws`, `kf_filter` (Kalman filter and RTS smoother),
`kf_stationary_cov`, `kf_forecast`, and `kf_fit_local_level`.

Probability vectors are columns and Markov transition matrices are row
stochastic. AR coefficients use recursion intercepts. Kalman initialization
is for the state before the first transition; missing measurements are omitted
from the observed-data likelihood. State intervals and forecasts condition on
fixed parameters. Each new help page includes a runnable standalone example.

Independent tests compare regime probabilities with exhaustive state-path
enumeration and Kalman results with dense joint-Gaussian conditioning, including
missing components, changing matrices and observation forecast variances.
Existing API and datasets are retained. No additional dependency or dataset
is introduced; the data provenance and licence register is unchanged.

# IdSS 0.2.2

Every exported function now has an executable help example. Added standalone
examples for the ICA objective and gradient, FEV-maximizing identification,
and variance decomposition. Examples are audited in separate fresh R sessions.
The ICA example also detected and corrects the Cayley-transform Jacobian used
by d.func.2.minimize; gradients are checked against numerical differences.

# IdSS 0.2.1

ARMAX regressors preserve their column names and accept explicit x_lags.
Coefficient tables no longer infer contemporaneous timing from column order.
Course compatibility tests cover the lag metadata. See estim.armax help.

The software remains MIT. JST remains CC BY-NC-SA 4.0, and USsentiment
retains noncommercial attribution conditions. See inst/DATA-SOURCES.md.
