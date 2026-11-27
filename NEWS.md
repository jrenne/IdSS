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
