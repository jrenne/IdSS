USmonetaryRR: frozen narrative-shock teaching data

Run build.R from the package root for an offline base-R reconstruction.
The CSV is checked by MD5; its SHA-256 and the original archive/source-file
SHA-256 are in provenance.json. Optional extract-from-upstream.R recreates
that CSV from the original author ZIP using haven, with no live download.
The 13 MB upstream ZIP is deliberately not included in the package/repository.

Source: https://github.com/ojorda/JEL-Code/tree/655696c1c576b7537c5a939d2c261f0a111ae663
Archive: LP_JEL_Replication.zip, Example6_JointInference/data_fred.dta.
The official author's repository applies CC0-1.0; its full notice and source
README are installed at data-sources/USmonetaryRR/. The Example 6 README
identifies the extended Romer-Romer measure and the public FRED inputs.

The five-column extraction is unchanged apart from date conversion, column
names and the documented 1969-03--2007-12 window. It excludes unrelated shock
series in the source file. No missing observation is filled, no interior month
is deleted. The macroeconomic inputs are archived Federal Reserve Board and
BLS observations, not a current-vintage reconstruction. Keep citations to
Jorda and Taylor (2025), Romer and Romer (2004), Coibion et al. (2017), the
original U.S. agencies and FRED. See ?USmonetaryRR for units and limitations.
