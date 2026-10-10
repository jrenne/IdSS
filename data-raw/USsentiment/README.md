# Rebuilding `USsentiment`

Run `Rscript data-raw/USsentiment/build_USsentiment.R` from the package root.
The script downloads current-vintage source files into a temporary directory,
constructs the six-column teaching dataset, checks its dimensions and missing
values, and writes `data/USsentiment.rda`.

The source files themselves are deliberately not vendored here. The generated
dataset consists of transformed economic observations from U.S. public-sector
sources. See `inst/DATA-SOURCES.md` for source terms and full attribution.

Because several inputs are revised over time, rebuilding in the future may not
be bit-for-bit identical to the packaged October 2026 vintage.
