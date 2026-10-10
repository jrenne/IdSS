# Data sources, transformations, and redistribution status

This file records the provenance and redistribution status of every dataset
shipped with IdSS. It is an audit of the package contents as of 10 October
2026; it is not legal advice. The MIT licence in `DESCRIPTION` applies to the
IdSS software, not automatically to third-party data.

Status labels used below:

- **CLEARED — attribution required**: the stated source terms permit
  redistribution, subject to the attribution shown below.
- **CLEARED — restricted licence**: redistribution is permitted only under the
  stated data licence.
- **HOLD — permission not documented**: the source is public, but no express
  redistribution licence was found. Do not include the dataset in a new public
  release until the rightsholder grants permission or an applicable licence is
  documented.
- **REMOVE/REPLACE**: the available evidence does not permit redistribution.

## Dataset register

| IdSS object | Contents and transformations | Original source and required attribution | Status |
|---|---|---|---|
| `BQ` | 159 quarterly U.S. observations, 1948 Q2–1987 Q4. Real-GNP log growth is computed from `GNPC96`, demeaned on either side of the historical break, and quarterly-average unemployment is linearly detrended. This is a teaching reconstruction, not the authors' original vintage dataset. | U.S. Bureau of Economic Analysis (`GNPC96`) and Bureau of Labor Statistics (`UNRATE`), retrieved through FRED. Cite both agencies and FRED. | **CLEARED — attribution required.** The inputs are public-domain U.S. government series. |
| `eurusd` | 6,770 non-missing daily observations, 1999-01-04–2025-12-31. FRED series `DEXUSEU` is renamed `USD_per_EUR`; no substantive transformation is applied. | Board of Governors of the Federal Reserve System, H.10 Foreign Exchange Rates, “U.S. Dollars to Euro Spot Exchange Rate [DEXUSEU],” retrieved from FRED, Federal Reserve Bank of St. Louis. | **CLEARED — attribution required.** FRED labels this series “Public Domain: Citation Requested.” |
| `international` | 279 quarterly observations, 1956 Q1–2025 Q3. Selected CPI, real-GDP/growth, and short-rate series were downloaded in a FRED workbook, renamed, and merged without substantive alteration. | OECD Main Economic Indicators via FRED. Cite the relevant OECD datasets and “retrieved from FRED, Federal Reserve Bank of St. Louis,” with an access date. | **CLEARED — attribution required.** OECD’s default data terms permit copying, adapting, distributing, and sharing for any purpose, unless series-specific or third-party restrictions apply. The archived workbook identifies these as OECD series. |
| `JST` | 2,718 annual country-year observations, 1870–2020; release 6 of the database is imported without substantive alteration. | Òscar Jordà, Moritz Schularick, and Alan M. Taylor (2017), “Macrofinancial History and the New Business Cycle Facts.” Additional citations are required for return and bank-balance-sheet variables; see the database site. | **CLEARED — restricted licence.** Creative Commons Attribution-NonCommercial-ShareAlike 4.0 (CC BY-NC-SA 4.0). Commercial use is not permitted, attribution is mandatory, and redistributed/adapted data must retain the same licence. |
| `USsentiment` | 168 quarterly observations, 1968 Q4–2010 Q3. Utilization-adjusted TFP growth is accumulated into a log index; nondurable-goods and services consumption is aggregated with a Törnqvist index and divided by population; real GDP and total hours are divided by population; the one-quarter-ahead Greenbook/Tealbook real-GDP-growth forecast is standardized. This is an independent current-vintage teaching reconstruction, not the authors’ replication dataset. | John Fernald’s quarterly TFP data (Federal Reserve Bank of San Francisco); BEA and BLS series `GDPC1`, `TOTLQ`, `DNDGRG3Q086SBEA`, `DSERRG3Q086SBEA`, `PCND`, `PCESV`, and `CNP16OV` via FRED; Greenbook/Tealbook projections created by Federal Reserve Board staff and transcribed/distributed by the Federal Reserve Bank of Philadelphia. Cite the agencies, data distributors, and Levchenko and Pandalai-Nayar (2020) as the motivating application. | **CLEARED — restricted licence.** BEA, BLS, and Board observations are U.S. public-sector data. The FRBSF website permits reproduction of its material without limitation in number when it is not for private gain and the Bank is credited. Accordingly, redistribution of this combined dataset is noncommercial and attribution is required. The package contains transformed observations, not the Philadelphia Fed workbook or the authors’ data. |
| `QuarterlyUS` | 314 quarterly observations, 1947 Q1–2025 Q2. Four FRED series are selected and renamed: `BOGZ1FL072052006Q`, `GDPC1`, `GDPDEF`, and `GDPPOT`. | Board of Governors of the Federal Reserve System, U.S. Bureau of Economic Analysis, and U.S. Congressional Budget Office, all retrieved through FRED. Cite each agency and FRED. | **CLEARED — attribution required.** These are U.S. government series; BEA and CBO explicitly identify their material as public domain, and FRED requests citation. |
| `US_CH_quarterly` | 117 quarterly observations, 1996 Q1–2025 Q1. Growth series derived by FRED from Swiss GDP (`CLVMNACSAB1GQCH`), Swiss and U.S. CPI (`CPALTT01CHQ659N`, `CPALTT01USQ659N`), and U.S. GDP (`GDPC1`). Columns are renamed. | Eurostat, OECD, and BEA via FRED. Cite each original producer and FRED and state that percentage-growth transformations were used. | **CLEARED — attribution required.** Eurostat authorises reuse with source acknowledgement; OECD permits distribution under its default data terms; BEA material is public domain. |
| `US3var` | 224 quarterly observations, 1959 Q2–2015 Q1. U.S. government series on GDP, potential GDP, unemployment, the natural rate, prices, the federal funds rate, and oil prices are merged; output and unemployment gaps, inflation, and commodity-price growth are computed. | BEA, BLS, CBO, Federal Reserve Board, and U.S. Energy Information Administration series archived from FRED. Cite the original agencies and FRED. | **CLEARED — attribution required**, provided the raw inputs remain the identified U.S. government series. |
| `USmonthly`, `USmonthlyExample`, `USmonthly_1990_to_2012` | Monthly U.S. activity, labour-market, price, monetary, reserve, interest-rate, and credit-risk series. Logs and first differences are computed where documented; three date windows are retained. `BRW` is the monthly Bu–Rogers–Wu monetary-policy shock and begins in 1994. | Board of Governors, BLS, and other U.S. government series retrieved through FRED; Bu, Rogers, and Wu shock series and Favara–Gilchrist–Lewis–Zakrajšek EBP downloaded from the Board of Governors. Cite the original series, FRED, and the corresponding research papers. | **CLEARED — attribution required.** The Board states that its website information is public domain unless otherwise indicated; the remaining inputs are U.S. government series. |
| `USyields` | 385 monthly observations, 1990–2022. Daily Treasury constant-maturity series `DGS1`, `DGS2`, `DGS3`, `DGS5`, `DGS7`, `DGS10`, `DGS20`, and `DGS30` are averaged by month. | Board of Governors of the Federal Reserve System, H.15 Selected Interest Rates, retrieved from FRED. | **CLEARED — attribution required.** FRED labels these series “Public Domain: Citation Requested.” |

## Source terms consulted

- FRED Services Terms of Use and FAQ:
  <https://fred.stlouisfed.org/legal>
- OECD Terms and Conditions, section “Data”:
  <https://www.oecd.org/en/about/terms-conditions.html>
- Eurostat copyright notice and free reuse of data:
  <https://ec.europa.eu/eurostat/help/copyright-notice>
- BEA public-domain statement:
  <https://www.bea.gov/help/faq/147>
- CBO copyright statement:
  <https://www.cbo.gov/about/privacy>
- Federal Reserve Board disclaimer and public-domain statement:
  <https://www.federalreserve.gov/disclaimer.htm>
- Federal Reserve Bank of San Francisco website policies (reproduction for
  noncommercial use with credit):
  <https://www.frbsf.org/site-policies/>
- Fernald quarterly TFP data:
  <https://www.frbsf.org/research-and-insights/data-and-indicators/total-factor-productivity-tfp/>
- Federal Reserve Board Greenbook/Tealbook data and Philadelphia Fed
  transcription documentation:
  <https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/greenbook>
- Bu--Rogers--Wu shock series:
  <https://www.federalreserve.gov/econres/feds/a-unified-measure-of-fed-monetary-policy-shocks.htm>
- Favara--Gilchrist--Lewis--Zakrajsek excess bond premium:
  <https://www.federalreserve.gov/econres/notes/feds-notes/updating-the-recession-risk-and-the-excess-bond-premium-20161006.html>
- JST database terms and required citations:
  <https://www.macrohistory.net/database/>

## Release rule

A public IdSS release should contain only datasets marked **CLEARED**. A
dataset marked **HOLD** may be restored to a public release after written
permission or an explicit licence is saved in `data-raw/<dataset>/` and this
register is updated. A dataset marked **REMOVE/REPLACE** must not be shipped.
The former `smi` dataset was replaced by the public-domain `eurusd` series. The
former `levpan`, Blanchard--Quah, and Ramey snapshots were replaced by
independent teaching reconstructions from redistributable government, Federal
Reserve, and attribution-compatible data. Future replacements must likewise
use a clearly licensed source, a simulated teaching dataset, or a user-side
download workflow whose terms permit that use.
