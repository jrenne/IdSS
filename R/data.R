#' Blanchard--Quah Quarterly Data
#'
#' Quarterly U.S. output-growth and detrended-unemployment series prepared for
#' the long-run-restrictions example of Blanchard and Quah.
#'
#' @format A data frame with 159 observations and three variables:
#' \describe{
#'   \item{Date}{Quarter date.}
#'   \item{Dgdp}{Real-GNP growth, in percent and adjusted for the sample break.}
#'   \item{unemp}{Linearly detrended unemployment rate.}
#' }
#' @source Blanchard--Quah replication data distributed through the
#'   Boston College Statistical Software Components archive,
#'   \url{https://econpapers.repec.org/software/bocbocode/rtz00017.htm}.
#' @references Blanchard, O. J. and Quah, D. (1989). The Dynamic Effects of
#'   Aggregate Demand and Supply Disturbances. *American Economic Review*,
#'   79(4), 655--673.
"BQ"

#' International Quarterly Macroeconomic Data
#'
#' Quarterly CPI and real-GDP series, together with selected short-term interest
#' rates, for Australia, Canada, Switzerland, the United Kingdom, Indonesia,
#' India, Japan, Korea, and the United States.
#'
#' @format A data frame with 279 observations and 35 variables. `date` is a
#'   `Date`. Country suffixes use ISO-like three-letter codes. Variables beginning
#'   with `CPI_`, `GDP_`, and `GDP1Y_` contain consumer prices, real GDP, and
#'   annual real-GDP growth; `STR_` variables contain short-term rates.
#' @source Course archive assembled from public macroeconomic series. See
#'   `data-raw/International/CPI_GDP.xlsx` in the source package for the archived
#'   input and series labels.
"international"

#' Jordà--Schularick--Taylor Macrohistory Database
#'
#' Long-run annual macroeconomic and financial data used for international
#' panel examples.
#'
#' @format A data frame with 2,718 country-year observations and 59 variables.
#'   Key identifiers are `year`, `country`, and `iso`. Remaining columns cover
#'   population, real activity, consumption, prices, trade, money, interest
#'   rates, housing, credit, public finance, crises, exchange-rate regimes, and
#'   total returns. Use `names(JST)` for the complete variable list.
#' @source Jordà--Schularick--Taylor Macrohistory Database,
#'   \url{https://www.macrohistory.net/database/}.
"JST"

#' Levchenko--Pandalai-Nayar Panel-VAR Replication Data
#'
#' Annual productivity, output, consumption, hours, and expected-duration data
#' used by the panel-VAR and FEVD examples.
#'
#' @format A data frame with 168 observations and six variables: `date`,
#'   `tfp_lev`, `lnconcap`, `lngdpcap`, `lnhrscap`, and `e12m`.
#' @source Replication archive supplied with the teaching materials; the original
#'   workbook is included as `data-raw/LevPan/levpan_replicationdata.xls`.
"levpan"

#' Quarterly U.S. Macroeconomic Levels
#'
#' A compact quarterly U.S. dataset containing the federal funds rate, real GDP,
#' the GDP deflator, and potential GDP.
#'
#' @format A data frame with 314 observations and five variables: `date`, `FFR`,
#'   `GDPC1`, `GDPDEF`, and `GDPPOT`.
#' @source Archived public U.S. macroeconomic series in
#'   `data-raw/QuarterlyUS/QVAR.xlsx`.
"QuarterlyUS"

#' Ramey Monetary-Policy and Macroeconomic Data
#'
#' Monthly U.S. macroeconomic variables and alternative monetary-policy shock
#' measures assembled for the local-projection and external-instrument examples.
#'
#' @format A data frame with 684 monthly observations and 42 variables. `DATES`
#'   is a `Date`. Macroeconomic variables include industrial production (`LIP`),
#'   unemployment (`UNEMP`), prices (`LCPI`, `LPCOM`), the federal funds rate
#'   (`FFR`), and the one-year Treasury rate (`GS1`). The remaining columns
#'   contain Romer--Romer, Bernanke--Kuttner, high-frequency, narrative, and
#'   factor-based shock series. Use `names(Ramey)` for the complete list.
#' @source Valerie A. Ramey's monetary-shock data archive,
#'   \url{https://econweb.ucsd.edu/~vramey/research.html#mon}.
"Ramey"

#' Quarterly U.S.--Swiss Macroeconomic Data
#'
#' Quarterly GDP and consumer-price series for Switzerland and the United States.
#'
#' @format A data frame with 117 observations and five variables: `date`,
#'   `CH_GDP`, `CH_CPI`, `US_CPI`, and `US_GDP`.
#' @source Archived public macroeconomic series in
#'   `data-raw/US_CH_quarterly/US_CH_quarterly.xlsx`.
"US_CH_quarterly"

#' Three-Variable U.S. VAR Dataset
#'
#' Quarterly U.S. series prepared for small structural-VAR examples.
#'
#' @format A data frame with 224 observations and six variables:
#' \describe{
#'   \item{Date}{Quarter date.}
#'   \item{y.gdp.gap}{Real-GDP gap, in percent.}
#'   \item{y.u.gap}{Unemployment gap, in percentage points.}
#'   \item{infl}{Annualized quarterly GDP-deflator inflation.}
#'   \item{r}{Federal funds rate.}
#'   \item{commo}{Commodity-price growth.}
#' }
#' @source Public U.S. macroeconomic series archived under
#'   `data-raw/US3var/` in the source package.
"US3var"

#' Monthly U.S. Macroeconomic and Monetary-Policy Data
#'
#' Monthly U.S. macroeconomic series merged with alternative monetary-policy
#' shock measures. Three sample variants are supplied for reproducible examples.
#'
#' @format `USmonthly` is a data frame with 670 observations and 49 variables;
#'   `USmonthlyExample` contains 366 observations; and
#'   `USmonthly_1990_to_2012` contains 266 observations. All have the same
#'   columns. `DATES` is a `Date`; core variables include `NBR`, `TTR`, `M1`,
#'   `LIP`, `UNEMP`, `LCPI`, `FFR`, `LPCOM`, `GS1`, and `EBP`. Other columns
#'   contain monetary-policy shock measures and factors. Use `names(USmonthly)`
#'   for the complete list.
#' @source Public macroeconomic series merged with Valerie A. Ramey's monetary
#'   policy data, \url{https://econweb.ucsd.edu/~vramey/research.html#mon}.
#' @name USmonthly
"USmonthly"

#' @rdname USmonthly
"USmonthlyExample"

#' @rdname USmonthly
"USmonthly_1990_to_2012"

#' Monthly U.S. Treasury Yields
#'
#' Monthly-average constant-maturity Treasury yields for eight maturities.
#'
#' @format A data frame with 385 observations and nine variables: `date`, `Y1`,
#'   `Y2`, `Y3`, `Y5`, `Y7`, `Y10`, `Y20`, and `Y30`. Yields are percentages and
#'   maturity suffixes are years.
#' @source Federal Reserve Economic Data (FRED), series `DGS1`, `DGS2`, `DGS3`,
#'   `DGS5`, `DGS7`, `DGS10`, `DGS20`, and `DGS30`,
#'   \url{https://fred.stlouisfed.org/}.
"USyields"

#' Swiss Market Index Daily Prices
#'
#' Daily open, high, low, close, adjusted-close, and volume observations for the
#' Swiss Market Index, used for the ARCH and GARCH examples.
#'
#' @format A data frame with 8,135 observations and seven variables:
#'   `Date`, `Open`, `High`, `Low`, `Close`, `Adj.Close`, and `Volume`. The sample
#'   runs from 9 November 1990 through 10 August 2022.
#' @source Archived market-data snapshot stored in `data-raw/SMI/SMI.csv` in the
#'   source package. Users should consult their data provider's terms before
#'   redistributing an updated snapshot.
"smi"
