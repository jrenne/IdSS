# Reading and maintaining IdSS

## Source map

| Files in `R/` | Purpose |
|---|---|
| `arma-garch.R` | ARMA/ARMAX estimation and simulation; ARCH/GARCH recursions |
| `var-simulation.R`, `shock-distributions.R` | VAR/VARMA simulation and independent shocks |
| `svar-recursive.R`, `svar-inference.R` | Recursive identification and sampling uncertainty |
| `svar-signs.R`, `svar-instruments.R`, `svar-fevmax.R` | Other identification methods |
| `hac-iv.R`, `local-projections.R` | HAC covariance, 2SLS, LP and LP-IV |
| `variance-decomposition.R` | FEVD calculations and plots |
| `ica.R`, `ica-distributions.R`, `ica-covariance.R` | ICA estimation, densities and covariance |
| `regime-switching.R` | Hamilton filtering, smoothing, simulation and forecasts |
| `kalman.R`, `latent-validation.R` | Kalman filtering, RTS smoothing and shared input checks |
| `data.R`, `data-licensing.R`, `IdSS-package.R` | Dataset and package documentation |

Public functions retain the names and arguments used by the companion courses.
Helpers stay beside their method; use descriptive names that do not resemble S3
methods of base generics. Keep matrix formulas close to their mathematical
notation, explain timing and normalization, and avoid disabled code blocks.
Historical experiments can be recovered from Git history.

## Style

Use two-space indentation, spaces around operators and after commas, one
statement per line, and short comments explaining the mathematics or a design
choice. Use `is.matrix()` and `is.list()` for type checks. Keep public argument
and result names stable, including their historical spelling. Write help in
roxygen comments; regenerate `man/` and `NAMESPACE` rather than editing them.
An optional development formatter is `styler::style_dir("R")`; styler is not a
package dependency.

## Checks

After `roxygen2::roxygenise()`, build and install the package into a temporary
library. Run `R CMD check` on the source archive. To run the tests against that
installed build from R:

```r
library(IdSS)
testthat::test_dir("tests/testthat", package = "IdSS",
  load_package = "installed", stop_on_failure = TRUE)
```

Run `Rscript --vanilla inst/validation/check-help-examples.R` using the same
library to audit every public example in a fresh session. Numerical changes
need independent benchmarks, and changes to the public API additionally require
rendering the companion courses. Keep data files and source-specific rights
notices separate from software refactoring; see `inst/DATA-SOURCES.md`.
