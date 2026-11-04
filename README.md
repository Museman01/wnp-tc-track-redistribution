# Western North Pacific tropical cyclone track redistribution

Analysis code for **Life-Cycle Redistribution of Western North Pacific Tropical
Cyclone Tracks and a Poleward Shift in Landfall Latitude: Cross-Agency Evidence
From 1966–2025**.

This study examines changes in tropical cyclone track composition and landfall
latitude over the western North Pacific. It compares the USA, Japan Meteorological
Agency (JMA), China Meteorological Administration (CMA), and Hong Kong Observatory
(HKO) records, separates genesis-distribution and conditional-track-field
contributions, and evaluates circulation, translation speed and upper-ocean
exposure along the tracks. Relative track changes are examined alongside absolute
passage, recurvature and landfall counts.

## Study scope

- The main comparison is **1966–1995 versus 1996–2025**. The relative
  track-composition domain is **100–180°E, 0–40°N**.
- Native-record comparisons retain each agency's selection criteria. The
  four-agency common sample contains **1,465 storms**.
- HKO is included throughout the 1966–2025 agency comparisons, with its official
  2025 post-analysis reports supplementing the IBTrACS record.
- The 1945-start comparison and 106-configuration robustness matrix retain their
  three-agency scope; HKO records begin in 1961.

## Main analyses

| Analysis | Main code |
| --- | --- |
| Cross-agency track fields, annual redistribution indices, common storms and sensitivities | `analysis/tracks/agency_comparison.py` |
| Genesis-distribution and conditional-track-field decomposition | `analysis/tracks/genesis_conditioned_decomposition.py`, `analysis/tracks/genesis_shapley.py` |
| Absolute passage, sustained recurvature and operational genesis | `analysis/tracks/track_baseline_diagnostics.py` |
| Landfall latitude, named-coast grouping and counting-rule sensitivities | `analysis/tracks/agency_comparison.py`, `analysis/landfall/hko_landfalls.py` |
| Circulation associations, deep/shallow steering and storm-sampling decomposition | `analysis/circulation/steering_layer_diagnostics.py`, `analysis/circulation/steering_sampling_decomposition.py`, `analysis/circulation/storm_time_sampling.py` |
| Trajectory-model calibration, baseline evaluation and finite sensitivities | `analysis/models/bam_primary_fit.py`, `analysis/models/bam_sensitivity.py` |
| Cross-agency translation-speed decomposition | `analysis/tracks/agency_comparison.py`, `analysis/tracks/speed_decomposition.py` |
| OHC300 and D26 exposure decomposition | `analysis/ocean/`, `analysis/ocean/d26_exposure.py` |
| Main and supplementary figures | `analysis/figures/main_figures.py`, `analysis/figures/si_four_agency_figures.py`, `analysis/figures/diagnostic_figures.py` |

Version 1.3 includes the four-agency comparison, absolute-count diagnostics,
trajectory-model and steering-layer sensitivities, the parallel D26 analysis, storm-matched six-hourly steering, daily circulation regimes, and count-based Table S7 reporting.
Statistical decompositions quantify contributions to observed differences;
circulation associations and model sensitivities are interpreted within their
respective observational and model assumptions.

## Code and data

The source code is organized by scientific topic, with shared utilities and
numerical tests. The accompanying science archive supplies the same source code,
scientific configurations, observations and regional input subsets, numerical
results, and seven main and twenty-one supplementary figures.

| Location | Contents |
| --- | --- |
| `analysis/tracks/` | Track composition, genesis, speed and pathway analyses |
| `analysis/landfall/` | Geometric crossings and coastal classification |
| `analysis/circulation/` | Circulation fields, steering and storm sampling |
| `analysis/ocean/` | OHC300 and D26 matching and exposure decomposition |
| `analysis/models/` | Trajectory-model calibration and sensitivity |
| `analysis/figures/` | Main and supplementary plotting |
| `analysis/shared/`, `core/` | Catalogue readers, statistics and shared utilities |
| `validation/tests/` | Numerical tests |
| `config/` | Scientific parameters |
| `data/` | Observations, processed tracks and regional fields; OHC300/D26 are in `data/derived/` |
| `results/` | All numerical outputs, grouped by scientific topic |
| `figures/` | Main and supplementary figures |
| `docs/` | Data sources, input methods and usage |

The included OHC300 and D26 regional fields cover 1966–2025. Three large ERA5
pressure-level and steering fields are obtained or regenerated separately before
recomputing fine-grid circulation, trajectory-model and storm-sampling analyses.
[docs/ERA5_INPUTS.md](docs/ERA5_INPUTS.md) specifies their source, variables, period, region,
grid, levels and processing methods; [docs/EXTERNAL_INPUTS.json](docs/EXTERNAL_INPUTS.json)
records the original file identifiers. Processed circulation arrays and results
are included in the science archive.

## Basic use

Use **Python 3.12**. Extract the science archive and run commands from its root;
it already contains the same Python sources as the source-only ZIP. Install the
dependencies listed in [requirements.txt](requirements.txt):

```bash
python -m pip install -r requirements.txt
```

The main settings are in `config/config.yaml`. Agency resampling settings are in
`config/agency_validation.yaml`; the ocean-exposure and pathway modules have
their own configuration files. Stage-specific seeds and resampling counts in the
source also form part of the analysis specification.

Use `python -m analysis --help` to list the analysis entry points.
All numerical outputs are organized under `results/`: `tracks/`, `agency_comparison/`,
`genesis/`, `landfall/`, `circulation/`, `models/`, `ocean/` and `pathways/`.
Diagnostic summaries are written to `results/diagnostics/`; pathway checks use
`results/pathways/diagnostics/`.

Examples using the distributed inputs:

```bash
# Cross-agency core statistics
python -B -m analysis tracks --execute core

# Landfall statistics from the included geometric crossing events
python -B -m analysis landfall --execute

# Main Figures 2 and 6
python -B -m analysis figures --execute 2 6

# Match tracks to the included regional OHC fields
python -B -m analysis ocean --phase prepare --use-prepared-ocean --execute --overwrite

# OHC decomposition and block-bootstrap intervals
python -B -m analysis ocean --phase analyze --execute --overwrite
```

These commands write tables or figures. Use a separate extraction to retain the
distributed results. See [docs/RUNNING.md](docs/RUNNING.md) for parameter locations and
input requirements. [requirements-download.txt](requirements-download.txt) lists
optional acquisition dependencies; Copernicus access requires a CDS account and
acceptance of the dataset terms.

## Data sources

- **IBTrACS v04r01:** [NOAA NCEI data record](https://doi.org/10.25921/82ty-9e16).
- **HKO 2025:** [Official post-analysis best-track data](https://data.weather.gov.hk/weatherAPI/hko_data/tc/HKO2025BST.csv).
- **ERA5 hourly pressure-level fields:** [C3S data record](https://doi.org/10.24381/cds.bd0915c6), sampled at 00/06/12/18 UTC.
- **ERA5 monthly pressure-level fields:** [C3S data record](https://doi.org/10.24381/cds.6860a573).
- **ORAS5 OHC300 and D26:** Copernicus Climate Change Service (2021).
  *ORAS5 global ocean reanalysis monthly data from 1958 to present* [Data set].
  ECMWF. [Dataset DOI](https://doi.org/10.24381/cds.67e8eeb7).

The science archive's `docs/SOURCES.md` also identifies the ONI/PDO indices, GSHHG
shorelines and relief-map inputs.

## Citation and licence

Please cite the software using [CITATION.cff](CITATION.cff), and cite the original
data products used in your analysis.

Huang, D., & Yue, Y. (2026). *Analysis code for post-genesis tropical-cyclone track
redistribution in the western North Pacific* (Version 1.3) [Software].
[Source repository](https://github.com/Museman01/wnp-tc-track-redistribution).

Software is released under the [MIT License](LICENSE). Project-generated
analytical tables and original figures use **CC BY 4.0**. Third-party
observations, reanalysis and geographic inputs retain their providers' terms;
these are described with the science archive and accompanying data notices.

## Storm-matched steering and daily circulation regimes

The additional circulation analysis uses the primary IBTrACS record with standard
LAT/LON coordinates and USA wind selection. Environmental steering is sampled
at six-hourly storm positions using a 5°–7° annulus and compared with the
same-year monthly mean at those positions. A 5° disk and a 3°–7° annulus provide
all-month comparisons. These are separate sampling geometries with their own
statistical families.

June–October daily winds are classified using a common 1966–2025 EOF basis and
four circulation regimes. Equal-year frequencies and persistence, fixed-region
winds, and the symmetric frequency/within-regime path decomposition support
Tables S13/S14 and Figures S20/S21. The shorter split is 1982–2003 versus
2004–2025. Three-year blocks, seed 202406, 9,999 permutations and 4,999 bootstrap
replicates specify the inference. Cluster-count sensitivity uses k=4–8,
17 principal components, 30 k-means initializations and seed 1.

From the science archive root:

```bash
python -B -m analysis.circulation.submonthly_analysis all --output submonthly_results --inference
python -B -m analysis.circulation.submonthly_products --output submonthly_products --refit-clusters
python -B -m analysis.figures.submonthly_steering . data/raw/GSHHG/GSHHS_shp/i/GSHHS_i_L1.shp submonthly_figures
python -B -m analysis.circulation.sampling_sensitivity --report-s7 --execute --output-csv table_s7.csv
```

`data/processed/submonthly_steering/` supplies matched wind vectors, standard
track points, daily regime labels, the circulation model and daily 2.5° winds.
`results/circulation/submonthly_steering/` contains statistical comparisons and
table data; `data/figure_inputs/submonthly_steering/` contains the plotting arrays.
The code ZIP contains Python sources; dependencies, scientific configurations and
data are provided by the science archive or accompanying repository files.

Statistics and plots can start from these processed inputs. Re-extracting storm
environmental winds from ERA5 requires the external six-hourly fields described
in [docs/SUBMONTHLY_STEERING.md](docs/SUBMONTHLY_STEERING.md). The distributed
wind statistics do not constitute four independent agency environmental analyses.
