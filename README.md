# Western North Pacific tropical cyclone track redistribution

Software version: **2.0.2**. Data archive: https://doi.org/10.5281/zenodo.24188084.

Analysis code and supporting data for *Life-Cycle Redistribution of Western
North Pacific Tropical Cyclone Tracks and a Poleward Shift in Landfall Latitude:
Cross-Agency Evidence From 1966–2025*.

The study compares USA, JMA, CMA and HKO tropical cyclone records, primarily
1966–1995 versus 1996–2025, over 100–180°E and 0–40°N. The analyses cover track
composition, genesis and conditional-track decomposition, landfall latitude,
circulation, recurvature, translation speed, trajectory-model sensitivities,
and exposure to OHC300 and the 26°C-isotherm depth.

## Contents

- `analysis/`: scientific calculations and main/supplementary figures.
- `core/`: shared paths, catalogue preparation and statistics.
- `config/`: scientific periods, grids, thresholds and other parameters.
- `data/`: observations, processed inputs and regional ocean/relief fields.
- `results/`: structured scientific estimates and diagnostics.
- `figures/`: six main and seventeen supplementary figures.

The code ZIP contains Python sources. The science archive also includes the
configurations, data, results, figures and documentation needed to use them.

## Run

Use Python 3.12 and run commands from the extracted science archive root:

```bash
python -m pip install -r requirements.txt
python -m analysis --help
python -B -m analysis tracks --execute core --output-root analysis_output/tracks
python -B -m analysis figures --execute 1 2 3 4 5 --output-dir analysis_output/figures
```

Choose new output locations to retain the distributed results.

| Documentation | Purpose |
| --- | --- |
| [RUNNING](docs/RUNNING.md) | Commands, parameter locations, scientific inputs and pre-landfall crossing criteria |
| [SOURCES](docs/SOURCES.md) | Data products, citations, licences and regional selections |
| [ERA5 inputs](docs/ERA5_INPUTS.md) | External pressure-level and steering inputs |
| [Circulation analysis](docs/SUBMONTHLY_STEERING.md) | Wind sampling and daily circulation regimes |
| [Additional diagnostics](docs/ADDITIONAL_DIAGNOSTICS.md) | Storm-level associations and height-gradient definitions |

Calculations and figures can start from the supplied processed inputs. Large
ERA5 pressure-level and six-hourly steering fields require separate acquisition
for meteorological preprocessing, environmental resampling and trajectory
integration. Their product selections and input formats are documented above.
Pre-landfall track truncation uses record-specific crossing criteria described
in RUNNING.

## Citation and licence

Cite the software using [CITATION.cff](CITATION.cff) and cite the original data
products used in your analysis. Source repository:
https://github.com/Museman01/wnp-tc-track-redistribution.

Software uses the [MIT License](LICENSE). Study-generated analytical tables
and original figures use CC BY 4.0. Third-party observations, reanalysis and
geographic inputs retain their providers’ terms.
