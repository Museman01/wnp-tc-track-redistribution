# Running the scientific code

## Software and parameters

Use Python 3.12 and install the accompanying dependencies:

```bash
python -m pip install -r requirements.txt
```

requirements-download.txt lists optional data-acquisition dependencies.
Copernicus access requires a CDS account and acceptance of the dataset terms.
The data archive supplies the data and configuration used by the source-only
ZIP. Run commands from the extracted data archive root.

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
python -B -m analysis figures --execute 1 2 3 4 5 --output-dir analysis_output/figures
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

The manuscript layouts for Figures 1–5 are drawn by
`analysis.figures.manuscript_figures`. The schematic and steering-layer
supplement use the separate commands in “Manuscript figures and index summaries”
below. Other plotting modules provide individual diagnostic components.

## Pre-landfall crossing criteria

In the three-agency pre-landfall sensitivity, tracks are truncated at the first
geometric shoreline crossing at any life-cycle stage. The crossing input is
`results/landfall/landfall_events_exact_full_lifecycle.csv`. With the distributed
native track records, the TV distances are 0.177707137275 (USA), 0.168856614191
(TOKYO/JMA) and 0.170746706482 (CMA), for 1966–1995 versus 1996–2025. These are
the values in
`data/reference/agency_comparison/04_track_density_sensitivity/path_definition_sensitivity.csv`.
Use this explicit input to calculate that definition:

```bash
python -B -m analysis.tracks.path_definition_validation --execute --events results/landfall/landfall_events_exact_full_lifecycle.csv --output-dir analysis_output/path_definitions_first_crossing
```

The default input, `results/landfall/landfall_events_exact.csv`, instead requires
the native report before a crossing to meet the agency's tropical-storm
threshold. Truncation at the first such qualified crossing yields TV distances
of 0.166910632903, 0.161948575851 and 0.163179752528 for USA, JMA and CMA,
respectively. The HKO pre-landfall sensitivity uses the first qualified
tropical-storm crossing. Under each criterion, storms without an eligible
crossing retain their full track. The native tropical-storm track-point
selection is the same for both truncation definitions; the crossing criterion
determines the cutoff time. Geometric landfall statistics elsewhere use the
qualified tropical-storm events.

## ERA5-dependent calculations

Obtain and prepare the pressure-level and steering inputs in docs/ERA5_INPUTS.md
before recomputing fine-grid circulation, trajectory-model or storm-sampling
analyses. Processed circulation arrays, tables and figures are included.
docs/EXTERNAL_INPUTS.json identifies the three original fields obtained or regenerated
separately; MANIFEST.csv identifies all distributed files.

## Manuscript figures and index summaries

Figures 1–5 use the structured track, circulation and landfall results. Figure 5
also reads `data/derived/etopo1_99-146E_4-46N_2min.nc`; selection and attribution
are in `docs/SOURCES.md`. The required scipy, xarray and geopandas packages are
listed in requirements.txt.

```bash
python -B -m analysis diagnostics --execute --only si --output-root analysis_output/diagnostic_figures
python -B -m analysis.figures.synthesis_figure --execute draw overlay analysis_output/Fig06_outlook_schematic
```

The primary/USA/JMA/CMA index summary and `index_hko_summary.csv` are stored in
`data/reference/environment_validation/supplemental_validation/results/`.
Use the `track_point` rows and the `late_minus_early_oos` and
`block_permutation_p` columns for the index differences and their p values.
The `tracks core` command with a new `--output-root` writes the HKO summary to
`reference_index/` within that output root and the other core tables to `core/`.


## Pathway sensitivity outputs

The pathway calculation writes structured membership, decomposition,
resampling and sensitivity tables under its selected output root. The fitted
classification is an analysis model; k/seed stability does not establish unique
physical pathway classes. Interpret the period contrast together with the
domain and density-denominator sensitivities. Statistical partitions and
conditional label-permutation tests do not identify dynamical causes.

Diagnostic counts, the trajectory-model baseline, steering-layer sensitivity
and D26 exposure are available through `python -m analysis diagnostics` with a
new `--output-root`. The D26 figure is S17. Figure 3 and Figure 1 contain the
circulation-geometry and four-agency map comparisons, respectively.
