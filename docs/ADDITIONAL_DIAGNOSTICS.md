# Additional track and circulation diagnostics

## Robustness matrix

`results/tracks/robustness_matrix_signed.csv` has 126 scientific configurations.
HKO contributes 20 configurations spanning track-point/storm-equal weighting,
three endpoint choices, and six leave-one-decade-out choices. The 2.5° domain is
100–180°E, 0–40°N. The signed projection uses the primary annual-field contrast
from 1966–1995 to 1996–2025. Three-year blocks do not cross gaps in retained years;
1,999 permutations use seed 202406. Undefined projections retain their stated
status. Recomputing HKO replaces the HKO records and preserves the non-HKO CSV
records and their stored precision.

## Storm-level circulation and recurvature

The anchor CSV has one row per storm: sid, season, anchor_time, anchor_method,
recurving, recurve_time, full_ts_reports and regime. Regime 0–3 corresponds to
regimes 1–4. Anchors use the full tropical-storm life history; eligible storms
supply in-domain tropical-storm points. Only June–October anchors enter the
1,155-storm sample. This retrospective anchor does not measure a forecast.

The regime–recurvature association uses a chi-square statistic and 9,999
within-season permutations of regime labels, with p=(b+1)/(B+1), seed 202406.
The logistic models compare a period indicator alone with period plus regime
dummies; the reported coefficient difference is adjusted minus unadjusted.
Intervals use 4,999 within-period nonoverlapping three-year block resamples;
9,999 three-year block permutations provide p values. Coefficients are log odds,
exponentiated coefficients are odds ratios, and standard errors cluster by year.
`q_bh_scalar_family` includes the four decomposition, four overall-period and
three logistic tests in the accompanying three CSV files. The supplied p values
for the other eight tests define the full BH family for the logistic output.

Persistence counts consecutive days of the anchor's regime through the anchor
itself. Counts start at 1 June when a run extends beyond the season boundary.
The weighted within-regime mean difference and its within-regime×season
permutation test use storm weights. A regime-adjusted logistic model uses log2
persistence; its odds ratio describes a doubling in elapsed persistence.

## Height gradients

`data/derived/circulation_monthly_coarse.npz` contains year, month, latitude,
longitude and z500 in metres. June–October annual mean fields have a 2.5° grid.
Centered differences use Earth radius 6,371 km and local spherical distances;
regional means use cosine-latitude weights. Regions are 110–180°E, 10–40°N
(domain); 120–160°E, 10–20°N (southern flank); 110–130°E, 15–30°N (western flank);
and 110–180°E, 25–40°N (northern side). The domain mean rise and spatial standard
deviation describe spatial structure; gradient units are metres per 1,000 km.
Period differences use 1966–1995 and 1996–2025. Moving three-year within-period
resamples provide intervals; nonoverlapping three-year block permutations provide
p values, with 4,999 and 9,999 draws respectively and seed 202406.

## Figures

Figure 3 uses the primary track-point contrast stored in
`results/tracks/track_point_contrast.npz`. Figure 7 separates the observed spatial
change, statistical partitions, interpretation, coastal/ocean expressions,
circulation association and untested processes. Its graphic retains
74%–80% as the single partition percentage in panel b. Plotting sources are under
`analysis/figures/`; all distributed numerical results retain their precision.

## Weighting and sample comparison

`results/genesis/field_comparison.csv` compares each record’s storm-equal shared-bin
contrast with its full annual-normalized track-point contrast. Spatial correlation
and projection onto the full contrast use the included spatial arrays. Run
`python -B -m analysis.tracks.field_comparison --output field_comparison.csv`.
