#' Open Reconstruction of the Blanchard--Quah Quarterly Example
#'
#' Quarterly U.S. real-GNP growth and detrended unemployment prepared from
#' public U.S. government series for the long-run-restrictions example of
#' Blanchard and Quah. This is a teaching reconstruction, not the authors'
#' original vintage dataset.
#'
#' @format A data frame with 159 observations and three variables:
#' \describe{
#'   \item{Date}{Quarter date.}
#'   \item{Dgdp}{Real-GNP growth, in percent and adjusted for the sample break.}
#'   \item{unemp}{Linearly detrended unemployment rate.}
#' }
#' @details **Redistribution status:** public-domain U.S. government series;
#'   attribution to the original agencies and FRED is requested. Quarterly real
#'   GNP growth is demeaned on either side of the historical sample break, and
#'   the quarterly-average unemployment rate is linearly detrended.
#' @source U.S. Bureau of Economic Analysis series `GNPC96` and U.S. Bureau of
#'   Labor Statistics series `UNRATE`, retrieved through FRED.
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
#' @details **Redistribution status:** permitted with attribution under the
#'   OECD data terms, subject to any series-specific or third-party restriction.
#'   The archived workbook identifies the included series as OECD series. See
#'   `DATA-SOURCES.md` in the installed package.
#' @source OECD Main Economic Indicators, retrieved through Federal Reserve
#'   Economic Data (FRED). The original FRED series identifiers are retained in
#'   `data-raw/International/CPI_GDP.xlsx` in the source package.
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
#' @details **Licence:** Creative Commons Attribution-NonCommercial-ShareAlike
#'   4.0 (CC BY-NC-SA 4.0). Use is non-commercial, attribution is mandatory,
#'   and shared copies or adaptations must use the same licence. The database
#'   site lists additional citations for return and bank-balance-sheet series.
#' @source Jordà--Schularick--Taylor Macrohistory Database,
#'   \url{https://www.macrohistory.net/database/}.
"JST"

#' Open U.S. Data for the Sentiment-Shock Example
#'
#' Quarterly U.S. productivity, output, consumption, hours, and expected-growth
#' data for a teaching reconstruction of the sentiment-shock example inspired
#' by Levchenko and Pandalai-Nayar (2020).
#'
#' @format A data frame with 168 quarterly observations and six variables:
#' \describe{
#'   \item{date}{Quarter date, from 1968 Q4 to 2010 Q3.}
#'   \item{tfp}{Log utilization-adjusted TFP index, accumulated from Fernald's
#'     annualized quarterly growth rate and normalized to zero in 1947 Q1.}
#'   \item{ln_consumption_pc}{Log Törnqvist index of real nondurable-goods and
#'     services consumption per capita, normalized to zero in 1947 Q1.}
#'   \item{ln_gdp_pc}{Log real GDP per member of the civilian
#'     noninstitutional population aged 16 and over, in chained dollars.}
#'   \item{ln_hours_pc}{Log annualized hours of wage and salary workers on
#'     nonfarm payrolls per member of that population.}
#'   \item{expectation}{Standardized one-quarter-ahead Greenbook/Tealbook
#'     forecast of annualized real-GDP growth.}
#' }
#' @details **Redistribution status:** cleared with attribution. The underlying
#'   BEA, BLS, and Federal Reserve Board observations are U.S. public-sector
#'   data. Fernald's TFP series may be reproduced for noncommercial use with
#'   source credit. The expectation series is a transformation of publicly
#'   released forecasts created by Federal Reserve Board staff. This is a
#'   current-vintage teaching reconstruction, not the authors' replication
#'   dataset. See `DATA-SOURCES.md` and the rebuild notes in `data-raw/`.
#' @source Fernald quarterly TFP data (Federal Reserve Bank of San Francisco);
#'   BEA and BLS series `GDPC1`, `TOTLQ`, `DNDGRG3Q086SBEA`,
#'   `DSERRG3Q086SBEA`, `PCND`, `PCESV`, and `CNP16OV`, retrieved through FRED;
#'   Federal Reserve Board Greenbook/Tealbook projections transcribed and
#'   distributed by the Federal Reserve Bank of Philadelphia.
#' @references Levchenko, A. A. and Pandalai-Nayar, N. (2020). TFP, News, and
#'   “Sentiments”: The International Transmission of Business Cycles. *Journal
#'   of the European Economic Association*, 18(1), 302--341.
"USsentiment"

#' Quarterly U.S. Macroeconomic Levels
#'
#' A compact quarterly U.S. dataset containing the federal funds rate, real GDP,
#' the GDP deflator, and potential GDP.
#'
#' @format A data frame with 314 observations and five variables: `date`, `FFR`,
#'   `GDPC1`, `GDPDEF`, and `GDPPOT`.
#' @details **Redistribution status:** public-domain U.S. government series;
#'   attribution to the original agencies and FRED is requested. `FFR` is FRED
#'   series `BOGZ1FL072052006Q`; the other series retain their FRED identifiers.
#' @source Board of Governors of the Federal Reserve System, U.S. Bureau of
#'   Economic Analysis, and U.S. Congressional Budget Office, retrieved through
#'   FRED. The archived workbook is `data-raw/QuarterlyUS/QVAR.xlsx`.
"QuarterlyUS"

#' Quarterly U.S.--Swiss Macroeconomic Data
#'
#' Quarterly GDP and consumer-price series for Switzerland and the United States.
#'
#' @format A data frame with 117 observations and five variables: `date`,
#'   `CH_GDP`, `CH_CPI`, `US_CPI`, and `US_GDP`.
#' @details **Redistribution status:** permitted with attribution under the
#'   Eurostat and OECD reuse terms; the BEA component is public domain. The four
#'   columns are FRED transformations of `CLVMNACSAB1GQCH`,
#'   `CPALTT01CHQ659N`, `CPALTT01USQ659N`, and `GDPC1`.
#' @source Eurostat, OECD, and U.S. Bureau of Economic Analysis, retrieved
#'   through FRED. Archived input: `data-raw/US_CH_quarterly/US_CH_quarterly.xlsx`.
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
#' @details **Redistribution status:** permitted with attribution. The source
#'   series are U.S. government data; IdSS merges them and computes output and
#'   unemployment gaps, inflation, and commodity-price growth.
#' @source U.S. Bureau of Economic Analysis, Bureau of Labor Statistics,
#'   Congressional Budget Office, Federal Reserve Board, and Energy Information
#'   Administration series retrieved through FRED and archived under
#'   `data-raw/US3var/`.
"US3var"

#' Monthly U.S. Macroeconomic and Monetary-Policy Data
#'
#' Monthly U.S. macroeconomic series and a redistributable monetary-policy shock.
#' Three sample variants are supplied for reproducible examples.
#'
#' @format `USmonthly` is a data frame with 670 observations and 16 variables;
#'   `USmonthlyExample` contains 366 observations; and
#'   `USmonthly_1990_to_2012` contains 266 observations. All have the same
#'   columns. `DATES` is a `Date`; core variables include `NBR`, `TTR`, `M1`,
#'   `LIP`, `UNEMP`, `LCPI`, `FFR`, `LPCOM`, `GS1`, and `EBP`. `BRW` is the
#'   monthly monetary-policy shock of Bu, Rogers, and Wu (2019); the remaining
#'   columns are first differences used by course examples.
#' @details **Redistribution status:** public-domain U.S. government and Federal
#'   Reserve Board series; attribution is requested. `BRW` begins in January
#'   1994, so earlier observations are missing. See `DATA-SOURCES.md` in the
#'   installed package.
#' @source Board of Governors of the Federal Reserve System, Bureau of Labor
#'   Statistics, and other U.S. government series retrieved through FRED; Bu,
#'   Rogers, and Wu monetary-policy shocks and the excess bond premium retrieved
#'   from the Board of Governors' Economic Research Data pages.
#' @references Bu, C., Rogers, J., and Wu, W. (2021). A Unified Measure of Fed
#'   Monetary Policy Shocks. *Journal of Monetary Economics*, 118, 331--349.
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
#' @details **Redistribution status:** FRED labels the underlying H.15 series
#'   “Public Domain: Citation Requested.” IdSS averages daily observations by
#'   calendar month.
#' @source Federal Reserve Economic Data (FRED), series `DGS1`, `DGS2`, `DGS3`,
#'   `DGS5`, `DGS7`, `DGS10`, `DGS20`, and `DGS30`,
#'   \url{https://fred.stlouisfed.org/}.
"USyields"

#' Daily U.S. Dollar–Euro Exchange Rate
#'
#' Daily noon buying rate in U.S. dollars per euro, used for the ARCH and GARCH
#' examples in the companion time-series course.
#'
#' @format A data frame with 6,770 observations and two variables:
#'   `Date`, the observation date, and `USD_per_EUR`, the number of U.S. dollars
#'   per euro. The sample runs from 4 January 1999 through 31 December 2025;
#'   dates with missing exchange-rate observations are omitted.
#' @details Daily percentage log returns can be computed as
#'   `100 * diff(log(eurusd$USD_per_EUR))`.
#'
#'   **Redistribution status:** FRED labels the underlying Board of Governors
#'   H.10 series “Public Domain: Citation Requested.”
#' @source Board of Governors of the Federal Reserve System (US), U.S. Dollars
#'   to Euro Spot Exchange Rate [DEXUSEU], retrieved from FRED, Federal Reserve
#'   Bank of St. Louis, \url{https://fred.stlouisfed.org/series/DEXUSEU}.
"eurusd"
