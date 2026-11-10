"""Association between the steering regime at the recurvature anchor and recurvature (June-October storms).

Input (read-only): storm_anchors_and_recurvature.csv from the ERA5 submonthly-steering analysis (stage 4),
one row per storm: season, regime label (0-3 = regimes 1-4 in Table S14), recurving (0/1).
Anchor: first native retained report at or north of 20N, otherwise the earliest westernmost report.
Test: chi-square statistic of the regime x recurvature table; null distribution from 9,999 permutations
of regime labels among storms within each season (keeps each season's regime and recurvature composition).
p = (b+1)/(B+1), seed 202406.
Usage: python -B regime_recurvature_association.py <anchors_csv> <out_json>
"""
import sys, json
import numpy as np, pandas as pd
src, out = sys.argv[1], sys.argv[2]
d = pd.read_csv(src)[['season', 'regime', 'recurving']].dropna()
d['regime'] = d.regime.astype(int); d['recurving'] = d.recurving.astype(int)
def chi2(reg, rec):
    t = np.zeros((4, 2))
    np.add.at(t, (reg, rec), 1)
    e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
    return float(((t - e) ** 2 / e).sum())
obs = chi2(d.regime.values, d.recurving.values)
rng = np.random.default_rng(202406)
groups = [np.flatnonzero(d.season.values == s) for s in np.unique(d.season.values)]
reg = d.regime.values.copy(); rec = d.recurving.values
B = 9999; b = 0
for _ in range(B):
    perm = reg.copy()
    for g in groups:
        perm[g] = reg[rng.permutation(g)]
    b += chi2(perm, rec) >= obs
p = (b + 1) / (B + 1)
tab = d.groupby('regime').recurving.agg(['count', 'sum', 'mean'])
res = dict(n_storms=int(len(d)), chi2=obs, permutation_p=p, n_permutations=B,
           by_regime={f'regime_{k+1}': dict(storms=int(r['count']), recurving=int(r['sum']), fraction=float(r['mean'])) for k, r in tab.iterrows()},
           early_storms=int((d.season <= 1995).sum()), late_storms=int((d.season >= 1996).sum()))
json.dump(res, open(out, 'w'), indent=2); print(json.dumps(res, indent=1))
