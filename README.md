# IdSS

`IdSS` provides the R procedures and datasets used by two companion courses:

- [Introduction to Time Series](https://jrenne.github.io/TimeSeries_Bookdown/);
- [The Identification of Dynamic Structural Shocks](https://jrenne.github.io/IdentifStructShocks/).

The package covers ARMA/ARMAX and ARCH/GARCH examples, VAR and VARMA
simulation, impulse responses, forecast-error variance decompositions, local
projections, and structural-shock identification through recursive orderings,
sign and zero restrictions, external instruments, FEVD maximization, and
non-Gaussian independent components. It also provides Markov-switching
simulation, filtering, smoothing, estimation and forecasting, and linear
Gaussian Kalman filtering, RTS smoothing and observation forecasts.

## Installation

Install the checked, tagged release:

```r
install.packages("remotes")
remotes::install_github("jrenne/IdSS", ref = "v0.3.2")
library(IdSS)
```

Every exported function and packaged dataset has an R help page. Start with:

```r
help(package = "IdSS")
?sim.arma
?svar.ordering
?USmonthly
?USmonetaryRR
?ms_filter
?fit_msar
?kf_filter
?kf_forecast
```

## Data provenance and licences

The package's MIT licence covers the IdSS software; it does not automatically
cover third-party data. The dataset-by-dataset register in
[`inst/DATA-SOURCES.md`](inst/DATA-SOURCES.md) records the original source,
transformations, required attribution, licence, and redistribution decision.

The former `BQ` replication snapshot and `Ramey`/`USmonthly` compilation have
been replaced by teaching datasets assembled from redistributable U.S.
government and Federal Reserve Board series. The uncleared SMI snapshot has
likewise been replaced by `eurusd`, a public-domain Federal Reserve H.10 series
retrieved through FRED. The uncleared `levpan` extract has been replaced by
`USsentiment`, an independently rebuilt U.S. dataset using public-sector data
and Fernald's noncommercial, attribution-required TFP series.

The `USmonetaryRR` dataset supplies the extended Romer--Romer instrument,
policy rate, unemployment and PCE inflation used by the LP-IV chapter. Its help
page documents all five variables, archived vintages, identification limits and
the official source repository's CC0 licence. Installed source evidence is at
`system.file("data-sources", "USmonetaryRR", package = "IdSS")`.

## Development checks

From the package directory:

```r
roxygen2::roxygenise()
library(IdSS)
testthat::test_dir("tests/testthat", package = "IdSS",
  load_package = "installed", stop_on_failure = TRUE)
devtools::check()
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the source map, coding conventions
and checks against an installed build.

Changes to the public API are additionally validated by rendering both companion
books against the candidate package version.

## Executable help examples

Every public function has a standalone example. After installing IdSS, run
`Rscript --vanilla inst/validation/check-help-examples.R` from the source
repository to execute each example in a separate fresh R session. CI runs
this audit as well as the ordinary R package examples and gradient checks.
