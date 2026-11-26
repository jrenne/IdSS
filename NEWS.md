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
