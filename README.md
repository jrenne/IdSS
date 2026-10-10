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

## Development checks

From the package directory:

```r
roxygen2::roxygenise()
devtools::test()
devtools::check()
```

Changes to the public API are additionally validated by rendering both companion
books against the candidate package version.
