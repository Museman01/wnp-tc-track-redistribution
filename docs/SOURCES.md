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
- ONI: https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso/oni/v5/
- PDO: https://www.ncei.noaa.gov/access/monitoring/pdo/
- GSHHG 2.3.7: https://doi.org/10.5281/zenodo.7007502
- Maps-For-Free relief: https://maps-for-free.com/

Original third-party terms apply to observations, reanalysis inputs, geographic
data and redistributed subsets. The project's MIT software license and CC BY
4.0 analytical-output license do not replace those terms. GSHHG notices are
retained with its input files.

## ETOPO1 relief input

Amante, C., & Eakins, B. W. (2009). *ETOPO1 1 Arc-Minute Global Relief Model:
Procedures, Data Sources and Analysis*. NOAA Technical Memorandum NESDIS
NGDC-24. https://doi.org/10.7289/V5C8276M

The grid-registered ice-surface altitude field was obtained on 2026-11-10 from
NOAA CoastWatch ERDDAP dataset `etopo180`, selecting 99–146°E and 4–46°N with
stride 2 on the native 1-arc-minute grid:
https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude%5B(4):2:(46)%5D%5B(99):2:(146)%5D

The supplied subset is `data/derived/etopo1_99-146E_4-46N_2min.nc`, with
1261 latitude and 1411 longitude nodes at 2-arc-minute spacing. `altitude` is
in metres, positive upward; horizontal datum WGS84 and vertical datum mean sea
level. The landfall figure uses shaded relief from this field and GSHHG L1
shorelines. NOAA's embedded terms permit free use and redistribution; the
original third-party terms remain applicable. A future download has not been
verified to match the supplied subset byte for byte.
