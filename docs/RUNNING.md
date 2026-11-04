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

## Storm-sampling permutation results

Table S7 uses p = (b + 1)/(B + 1), where b is the number of permutation
statistics at least as extreme as the observed statistic and B = 9,999.
USA, JMA and CMA counts are supplied in
`data/reference/storm_sampling/permutation_test.json`; HKO counts are in
`results/agency_comparison/table_s7_permutation_four_agency.csv`.
The resulting table is `results/circulation/storm_sampling/table_s7_reporting.csv`.

```bash
python -B -m analysis.circulation.sampling_sensitivity --report-s7 --execute --output-csv s7_reporting.csv
```

The command uses the supplied exceedance counts. Values below 0.001 are
reported as <0.001; other displayed p values are rounded to three decimal places.

## Storm-matched steering and regimes

See [SUBMONTHLY_STEERING.md](SUBMONTHLY_STEERING.md) for processed inputs, external ERA5 selection and scientific parameters. The README gives commands for the statistics, Tables S13/S14 and Figures S20/S21.
