# Running the scientific code

## Software and parameters

Use Python 3.12 and install the accompanying dependencies:

```bash
python -m pip install -r requirements.txt
```

requirements-download.txt lists optional data-acquisition dependencies.
Copernicus access requires a CDS account and acceptance of the dataset terms.
The science archive supplies the data and configuration used by the source-only
ZIP. Run commands from the extracted science archive root.

- config/config.yaml: domains, periods, intensity thresholds, season and grids.
- config/agency_validation.yaml: agency resampling and exact-landfall settings.
- config/ohc.json: ocean matching, year assignment,
  land treatment, weighting, product boundary and block bootstrap.
- config/pathways.json: pathway classification and decomposition.

Stage-specific random seeds and resampling counts in the source are also part
of the analysis specification.

## Selected calculations

Choose new output locations to retain the distributed results. These examples
use included observations, processed events, numerical fields or regional ocean
inputs:

```bash
python -B -m analysis tracks --execute core --output-root analysis_output/agency
python -B -m analysis landfall --execute --output-root analysis_output/agency
python -B -m analysis figures --execute 2 6 --output-dir analysis_output/figures
python -B -m analysis.ocean.ohc_summary --execute --output-dir analysis_output/ohc_summary
```

The landfall calculation starts from the distributed geometric crossing events.
To recalculate native-stage landfall statistics in a separate output directory:

```bash
python -B -m analysis.landfall.native_landfall_validation --execute --reuse-existing --output-dir analysis_output/native_landfall
```

To refine new crossings instead, supply `--candidate-events` with a CSV of
candidate crossings and omit `--reuse-existing`.
The OHC summary uses the included matched-point and decomposition results. OHC
preparation and decomposition use the paths and settings in `config/ohc.json`.

ORAS5 uses the consolidated product through 2014 and the operational product
from 2015. The analysis retains these fields without automatic bias correction.

When the original D26 archives and their request records are available, regional
fields or copies with native timestamps can be written to new directories:

```bash
python -B -m analysis.ocean.d26_preparation --execute --output-dir analysis_output/d26_fields
python -B -m analysis.ocean.d26_time_coordinates --execute --output-dir analysis_output/d26_times
```

The second command reads the existing regional fields and writes copies; it does
not modify the original fields.


HKO figure-series tables and the study-region map can also be written separately:

```bash
python -B -m analysis.figures.hko_figure_series basic --execute --output-dir analysis_output/hko_series
python -B -m analysis.figures.study_region_figure --execute --output analysis_output/study_region
```

The figure command writes PDF and PNG files using the requested output stem.


## Standalone decomposition and sensitivity calculations

These entry points use the distributed track records, exact landfall events,
annual-composition arrays and annual core metrics. Each calculation requires a
new output directory. They retain the three-agency sample definitions of the
corresponding analyses; four-agency comparisons use `analysis tracks`.

```bash
python -B -m analysis.tracks.agency_composition_validation --execute --output-dir analysis_output/three_agency_core
python -B -m analysis.tracks.common_storm_validation --execute --output-dir analysis_output/common_storms
python -B -m analysis.tracks.climate_adjustment_validation --execute --output-dir analysis_output/climate_adjustment
python -B -m analysis.tracks.cutpoint_validation --execute --output-dir analysis_output/cutpoints
python -B -m analysis.tracks.path_definition_validation --execute --output-dir analysis_output/path_definitions
python -B -m analysis.landfall.coast_grouping_validation --execute --output-dir analysis_output/coast_grouping
python -B -m analysis genesis PRIMARY --execute --output-dir analysis_output/genesis_primary
python -B -m analysis.tracks.genesis_shapley --execute --output-dir analysis_output/genesis_records
python -B -m analysis.tracks.spatial_redistribution --execute --output-dir analysis_output/spatial_redistribution
```

Landfall-based sensitivity calculations accept `--events`. Climate adjustment
and cutpoint calculations accept `--composition`; cutpoint calculations also
accept `--core-annual`. These options select scientific input files and leave
the configured periods and statistical parameters unchanged.
Cutpoint and spatial-redistribution calculations write numerical results by
default; add `--plot` to request their diagnostic figures.

For uniform lettering of supplementary diagnostic figures S14–S19, select the
figure number and a new output directory, for example:

```bash
python -B -m analysis.figures.uniform_typography S14 --execute --output-dir analysis_output/figure_s14
```

## Pre-landfall reference-table limitation

The retained three-agency table
`data/reference/agency_comparison/04_track_density_sensitivity/path_definition_sensitivity.csv`
contains pre-landfall TVDs of 0.177707137275 (USA), 0.168856614191 (TOKYO/JMA)
and 0.170746706482 (CMA). The corresponding event-input provenance is not established.
With the distributed `results/landfall/landfall_events_exact.csv` and track records,
the pre-landfall definition yields TVDs of 0.166910632903, 0.161948575851 and
0.163179752528, respectively, for 1966–1995 versus 1996–2025.
The retained reference values should therefore not be used as expected outputs
for these supplied events. These point estimates do not establish significance
levels; inference depends on the specified resampling procedure.

## ERA5-dependent calculations

Obtain and prepare the pressure-level and steering inputs in docs/ERA5_INPUTS.md
before recomputing fine-grid circulation, trajectory-model or storm-sampling
analyses. Processed circulation arrays, tables and figures are included.
docs/EXTERNAL_INPUTS.json identifies the three original fields obtained or regenerated
separately; MANIFEST.csv identifies all distributed files.
