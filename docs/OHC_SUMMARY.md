# OHC exposure summary tables

`analysis/ocean/ohc_summary.py` reads the matched-point, decomposition and product-boundary results from `results/ocean/ohc300/` and writes main estimates, configuration counts, factor sensitivities and spatial/temporal summaries. Estimates remain in J m-2 unless a column states another unit.

```bash
python -B -m analysis.ocean.ohc_summary --execute --output-dir analysis_output/ohc_summary
```
