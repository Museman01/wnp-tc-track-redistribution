# Scientific data sources

- IBTrACS Western Pacific v04r01: https://doi.org/10.25921/82ty-9e16
- HKO 2025 post-analysis reports:
  https://data.weather.gov.hk/weatherAPI/hko_data/tc/HKO2025BST.csv
- ERA5 monthly pressure-level fields: https://doi.org/10.24381/cds.6860a573
  Selection, normalization and derived steering fields: docs/ERA5_INPUTS.md;
  original file identifiers: docs/EXTERNAL_INPUTS.json.
- ORAS5 OHC300 and D26 monthly ocean fields: Copernicus Climate Change Service
  (2021). *ORAS5 global ocean reanalysis monthly data from 1958 to present*
  [Data set]. ECMWF. https://doi.org/10.24381/cds.67e8eeb7
  Access: https://cds.climate.copernicus.eu/datasets/reanalysis-oras5
- ONI: https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
- PDO: https://www.ncei.noaa.gov/access/monitoring/pdo/
- GSHHG 2.3.7: https://doi.org/10.5281/zenodo.7007502
- Maps-For-Free relief: https://maps-for-free.com/

Original third-party terms apply to observations, reanalysis inputs, geographic
data and redistributed subsets. The project's MIT software license and CC BY
4.0 analytical-output license do not replace those terms. GSHHG notices are
retained with its input files.

Table S7 reporting data are in
`results/circulation/storm_sampling/table_s7_reporting.csv`. They combine the
retained three-agency exceedance counts in
`data/reference/storm_sampling/permutation_test.json` with the HKO counts in
`results/agency_comparison/table_s7_permutation_four_agency.csv`.
The corrected three-agency values are also supplied in
`results/circulation/storm_sampling/permutation_test.json`.
The reporting formula and the one-count differences from the vectorized
comparison are described in `docs/RUNNING.md`.
