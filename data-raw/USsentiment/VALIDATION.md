# Validation against the former `levpan` extract

The October 2026 reconstruction was compared with the former 168-quarter
`levpan` extract over 1968 Q4--2010 Q3. The old file was used only for this
validation and was not used to construct or normalize `USsentiment`.

For the four macroeconomic series, correlations with the old extract were:

| Series | Levels | First differences |
|---|---:|---:|
| Utilization-adjusted TFP | 0.9997 | 0.9866 |
| Real GDP per capita | 0.9999 | 0.9966 |
| Real nondurable-and-services consumption per capita | 0.9997 | 0.8891 |
| Hours per capita | 0.9990 | 0.9756 |

The expectation variable is intentionally different: a standardized
one-quarter-ahead Greenbook/Tealbook real-GDP-growth forecast replaces the
standardized Michigan E12 sentiment index in the old extract.

The course's five-variable VAR(2), FEVD-max identification, and 40-quarter
impulse responses were then re-estimated with the same random seed and 100
bootstrap replications. For the identified sentiment shock, correlations
between the old and reconstructed median impulse-response paths were 0.813
(TFP), 0.960 (GDP), 0.960 (expectation/sentiment measure), 0.934 (consumption),
and 0.958 (hours). At horizon 40, the sentiment shock accounts for 12.7 percent
of the bootstrap-average GDP forecast-error variance with the old data and
10.0 percent with the reconstruction; the combined news-and-sentiment share is
72.7 and 69.6 percent, respectively. The qualitative result and relative
persistence are preserved. This confirms that `USsentiment` is suitable as a
close teaching reconstruction, but it should not be presented as an exact
replication of the authors' data or estimates.
