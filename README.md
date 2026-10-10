# IdSS

`IdSS` provides the R procedures and datasets used by two companion courses:

- [Introduction to Time Series](https://jrenne.github.io/TimeSeries_Bookdown/);
- [The Identification of Dynamic Structural Shocks](https://jrenne.github.io/IdentifStructShocks/).

The package covers ARMA/ARMAX and ARCH/GARCH examples, VAR and VARMA
simulation, impulse responses, forecast-error variance decompositions, local
projections, and structural-shock identification through recursive orderings,
sign and zero restrictions, external instruments, FEVD maximization, and
non-Gaussian independent components.

## Installation

```r
install.packages("remotes")
remotes::install_github("jrenne/IdSS")
library(IdSS)
```

Every exported function and packaged dataset has an R help page. Start with:

```r
help(package = "IdSS")
?sim.arma
?svar.ordering
?USmonthly
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

## Development checks

From the package directory:

```r
roxygen2::roxygenise()
devtools::test()
devtools::check()
```

Changes to the public API are additionally validated by rendering both companion
books against the candidate package version.
