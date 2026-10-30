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
- config/agency_validation.yaml: agency resampling settings.
- config/ohc.json: ocean matching, year assignment,
  land treatment, weighting, product boundary and block bootstrap.
- config/pathways.json: pathway classification and decomposition.

Stage-specific random seeds and resampling counts in the source are also part
of the analysis specification.

## Selected calculations

Run in a separate extraction to retain the distributed results. These examples
use included observations, processed events, numerical fields or regional ocean
inputs:

```bash
python -B -m analysis tracks --execute core
python -B -m analysis landfall --execute
python -B -m analysis figures --execute 2 6
python -B -m analysis ocean --phase prepare --use-prepared-ocean --execute --overwrite
python -B -m analysis ocean --phase analyze --execute --overwrite
```

The landfall calculation starts from the distributed geometric crossing events.
The OHC preparation command extracts tracks and matches them to the included
regional ocean fields. The analysis command runs the decomposition configurations
and bootstrap replicates specified in the OHC configuration.

ORAS5 uses the consolidated product through 2014 and the operational product
from 2015. The analysis retains these fields without automatic bias correction.

## ERA5-dependent calculations

Obtain and prepare the pressure-level and steering inputs in docs/ERA5_INPUTS.md
before recomputing fine-grid circulation, trajectory-model or storm-sampling
analyses. Processed circulation arrays, tables and figures are included.
docs/EXTERNAL_INPUTS.json identifies the three original fields obtained or regenerated
separately; MANIFEST.csv identifies all distributed files.
## Table S7 reporting

```bash
python -B -m analysis.circulation.sampling_sensitivity --report-s7 --execute --output-csv s7_reporting.csv
```

The command reads included records and writes a new CSV. It does not require
ERA5 fields or generate a new permutation sample. Existing output files are
preserved unless `--overwrite` is supplied.

For USA, JMA and CMA, integer exceedance counts are recovered from the original
`data/reference/storm_sampling/permutation_test.json` values, which store b/9999.
Their corrected values are checked against
`results/circulation/storm_sampling/permutation_test.json`. HKO counts come from
`results/agency_comparison/table_s7_permutation_four_agency.csv`.
All reporting p values use (b+1)/10000. The CSV retains the counts, four-decimal
p values and displayed values: values below 0.001 are shown as <0.001; other
values use decimal half-up rounding to three places.

The vectorized diagnostic has one-exceedance differences in six CMA zonal terms.
Table S7 retains the original three-agency counts throughout. Thus its CMA zonal
fixed-margin interannual p is 0.4704 (display 0.470), while the vectorized
calculation gives 0.4705. Both records remain available for comparison.
