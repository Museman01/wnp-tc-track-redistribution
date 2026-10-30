# OHC exposure summary tables

analysis/ocean/ohc_summary.py reads the matched-point, decomposition and
product-boundary results from results/ocean/ohc300/ and writes main estimates,
configuration counts, factor sensitivities and spatial/temporal summaries
in that same result directory. Estimates remain in J m-2 unless a table column
explicitly states another unit.

Run from the archive root:

```bash
python -B -m analysis.ocean.ohc_summary --execute --overwrite
```
