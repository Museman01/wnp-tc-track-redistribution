# ERA5 inputs and steering methods

## Source and selection

Copernicus Climate Change Service, ERA5 monthly averaged data on pressure levels
from 1940 to present: https://doi.org/10.24381/cds.6860a573 .
Access: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels-monthly-means .
Select monthly means of daily means (monthly_averaged_reanalysis), deterministic
final reanalysis (expver 0001, member 0 where those coordinates are supplied).
A CDS account and acceptance of the dataset licence are required.

| Selection | Value |
| --- | --- |
| Period | January 1940–December 2025, all months; 1,032 monthly records |
| Region | 80–180°E, 0–65°N |
| Grid | 0.25°; latitude 65 down to 0 (261), longitude 80 up to 180 (401) |
| Wind | u_component_of_wind and v_component_of_wind, m s-1 |
| Wind levels | 850, 700, 600, 500, 400, 300, 200 hPa |
| Geopotential | geopotential at 500 hPa, m2 s-2 |
| Relative vorticity | vorticity at 850 hPa, s-1 |

## Pressure-level preparation

The analysis reads data/interim/era5_wnp_dynamic_plev.nc, with variables
u, v, z and vo on (time, level, latitude, longitude), shape
(1032, 7, 261, 401), stored as float32. Where necessary rename valid_time
to time and pressure_level or isobaricInhPa to level. Combine the downloaded
months in chronological order, with one record per month and the coordinate
order given above. No spatial interpolation, temporal smoothing or unit
conversion is applied. Retain z as geopotential, rather than geopotential height.

Use level order [850, 700, 600, 500, 400, 300, 200] hPa. Wind has all seven
levels; z is populated only at 500 hPa and vo only at 850 hPa. Other levels
of z and vo are missing, not zero. Use month-start timestamps at 00:00;
the original file uses seconds since 1970-01-01 and proleptic_gregorian.
The requested wind slices, z500 and vo850 must be finite on the selected grid.

## Steering derivation

From the extracted archive root, after preparing the pressure-level input:

```bash
python -B -m analysis.circulation.steering_flow --output-dir analysis_output/era5
python -B -c "from analysis.circulation.steering_layer_diagnostics import derive; derive(steering_path='analysis_output/era5/steering.nc', output_dir='analysis_output/era5')"
```

The deep-layer wind is the trapezoidal pressure integral at 300, 400, 500,
600, 700 and 850 hPa divided by 550 hPa, for each wind component.
Shear is sqrt((u200-u850)^2 + (v200-v850)^2).
data/interim/steering.nc contains u_steer, v_steer and shear for 1940–2025,
shape (1032, 261, 401), in m s-1. The first command also produces
data/processed/steering_annual.csv.

The shallow-layer wind uses [500, 600, 700, 850] hPa and trapezoidal weights
[50, 100, 125, 75]/350. data/derived/steering_850_500.nc contains
u_steer and v_steer for 1966–2025, shape (720, 261, 401), in m s-1.
The second command also generates data/derived/circulation_monthly_coarse.npz.
These commands write the regenerated fields to `analysis_output/era5/`. Analyses
that read the canonical `data/interim/` or `data/derived/` paths require the
corresponding inputs at those locations.

The three pressure-level/steering files are obtained or regenerated separately
for fine-grid ERA5 analyses. Regional OHC300/D26 inputs and processed circulation
results remain included. docs/EXTERNAL_INPUTS.json records the original field sizes
and hashes; a new provider export can differ in encoding or metadata and may
therefore have a different file hash.
