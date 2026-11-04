# Submonthly steering inputs and calculations

## Data and sampling

The meteorological source is Copernicus Climate Change Service (2018), *ERA5 hourly
data on pressure levels from 1940 to present*, https://doi.org/10.24381/cds.bd0915c6.
The catalogue is https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels.
Retrieval requires a CDS account and acceptance of the dataset licence. Select
reanalysis u and v at 850, 700, 600, 500, 400 and 300 hPa, and geopotential at
500 hPa; 1966–2025, all months and valid days, 00/06/12/18 UTC; the box
60°N–10°S, 90–190°E; a regular 1° grid; NetCDF format. The 850–500-hPa comparison
uses 850, 700, 600 and 500 hPa. Pressure-weighted steering uses trapezoidal
integration over pressure divided by layer thickness.

An annual input for environmental resampling has coordinates time, latitude
(60 to −10 in descending 1° steps), and longitude (90 to 190 in ascending 1°
steps). Fields u_deep, v_deep, u_shallow, v_shallow and u500 are in m s⁻¹;
z500 is geopotential in m² s⁻². Time contains every four-times-daily UTC sample
in that calendar year. An external-input location is
`data/external/era5/six_hourly/<year>.nc`. These large annual fields are not
included. `geometry_weights`, `sampled_month_mean` and `extract_point` in
`analysis/circulation/submonthly_statistics.py` define the spatial masks,
cosine-latitude weighting, complete-geometry requirement and monthly comparisons.
Monthly wind here means the average of the same year's four-times-daily samples.

The matched-wind table uses standard IBTrACS coordinates, USA wind selection,
TS-stage points at or above 34 kt, and successive reports no more than 12 h
apart. Its 29,557 displacement segments retain the 180°E edge where the expanded
ERA5 grid supports sampling. The 5°–7° annulus is evaluated in June–October;
the 5° disk and 3°–7° annulus are all-month comparisons. The geometry prefixes
in the table are annulus_5_7, disk_5 and annulus_3_7. Wind components are m s⁻¹;
motion minus environmental wind defines residual motion.

## Daily circulation and path composition

Daily winds average the 00/06/12/18 UTC samples. Spherical cell overlap maps the
deep-layer winds to 2.5° cells spanning 0–45°N, 100–180°E. Each annual NPZ in
`data/processed/submonthly_steering/daily_winds/` contains 153 June–October days:
field (153,2,18,32), ridge (153), favorable (153), and fixed (153,3).
The ridge is the 500-hPa easterly-to-westerly transition over 110–150°E;
favorable is the sign of zonal steering at 25°N averaged over 125–145°E;
fixed contains cosine-weighted u, v and the positive-u indicator for 10–30°N,
110–160°E. These arrays and the matched vectors are processed starting inputs.

EOF anomalies remove the pooled 60-year calendar-day climatology and apply
sqrt(cos(latitude)) weights. The saved circulation model retains 17 components
and four regimes; labels are ordered by the first-PC cluster-centroid coordinate.
Persistence is the mean run length within each June–October season, with runs
censored at the season boundaries. Frequency and persistence comparisons give
each year equal weight. Path compositions use annual-normalized June–October
points on the 2.5° 0–40°N, 100–180°E grid and match the daily regime labels.

The symmetric decomposition separates changes in regime frequency from changes
in conditional path composition within regimes. Projection percentages use the
observed spatial contrast as reference. Three-year blocks, seed 202406,
9,999 permutations and 4,999 bootstrap replicates specify inference. For the
44-year comparison, each 22-year period contains seven triples and one singleton;
permutations exchange blocks within their lengths. Wind families and regime/
physical-index families are identified in the corresponding CSV columns.

The included model, labels, daily fields, matched vectors, standard track tables,
comparison tables and figure arrays define the processed-data calculations.
Fresh CDS retrievals can reflect later data updates and are not assumed to be
byte-identical to the inputs used here. Re-extracting meteorological fields or
reintegrating trajectories requires those external fields; trajectories driven
by storm-center winds are not part of the results presented here.
