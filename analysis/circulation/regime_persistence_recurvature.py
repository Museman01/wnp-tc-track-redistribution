"""Association between steering-regime persistence at the recurvature anchor and recurvature (June-October storms).

Inputs (read-only):
  anchors CSV  - storm_anchors_and_recurvature.csv (season, anchor_time, regime 0-3 = regimes 1-4, recurving 0/1)
  regimes CSV  - stage3_daily_regimes.csv (date, regime) for all 9,180 June-October days, 1966-2025
Persistence: number of consecutive days, up to and including the anchor day, on which the regime equals the anchor-day
regime. Runs do not cross a season boundary, so runs that start on 1 June are left-censored (counted from 1 June).
Only days up to the anchor day enter; no later information is used.
Statistic: within-regime difference in mean persistence, recurving minus nonrecurving storms, averaged over the four
regimes with weights equal to the number of storms in each regime. Null distribution: 9,999 permutations of persistence
values among storms of the same regime and season (keeps the regime-recurvature association and each season's regime
composition). Two-sided p = (b+1)/(B+1), seed 202406. A storm-weighted logistic model of recurvature on regime and
log(persistence) gives an odds ratio per doubling of persistence; its interval comes from 4,999 resamples of
nonoverlapping three-year blocks of seasons.
Usage: python -B regime_persistence_recurvature.py <anchors_csv> <regimes_csv> <out_json>
"""
import sys, json
import numpy as np, pandas as pd
anc, reg, out = sys.argv[1:4]
a = pd.read_csv(anc); r = pd.read_csv(reg)
r['date'] = pd.to_datetime(r['date']).dt.date
r = r.sort_values('date').reset_index(drop=True)
lab = dict(zip(r.date, r.regime.astype(int)))
# elapsed run length ending at each day
run = {}; prev_d = None; prev_l = None; cnt = 0
for d, l in zip(r.date, r.regime.astype(int)):
    if prev_d is not None and (pd.Timestamp(d) - pd.Timestamp(prev_d)).days == 1 and l == prev_l:
        cnt += 1
    else:
        cnt = 1
    run[d] = cnt; prev_d, prev_l = d, l
a['adate'] = pd.to_datetime(a.anchor_time).dt.tz_localize(None).dt.date if pd.to_datetime(a.anchor_time).dt.tz is not None else pd.to_datetime(a.anchor_time).dt.date
a = a.dropna(subset=['regime'])
a['regime'] = a.regime.astype(int); a['recurving'] = a.recurving.astype(int)
assert all(lab.get(d) == g for d, g in zip(a.adate, a.regime)), 'anchor regime differs from daily label'
a['pers'] = [run[d] for d in a.adate]
def stat(p, rec, rg):
    tot = 0.0; n = 0
    for k in range(4):
        m = rg == k
        if rec[m].min() == rec[m].max(): continue
        tot += m.sum() * (p[m & (rec == 1)].mean() - p[m & (rec == 0)].mean()); n += m.sum()
    return tot / n
p = a.pers.values.astype(float); rec = a.recurving.values; rg = a.regime.values; ss = a.season.values
obs = stat(p, rec, rg)
rng = np.random.default_rng(202406)
strata = [np.flatnonzero((rg == k) & (ss == s)) for k in range(4) for s in np.unique(ss)]
strata = [g for g in strata if len(g) > 1]
B = 9999; b = 0
for _ in range(B):
    q = p.copy()
    for g in strata: q[g] = p[rng.permutation(g)]
    b += abs(stat(q, rec, rg)) >= abs(obs)
pval = (b + 1) / (B + 1)
# logistic: recurving ~ regime dummies + log2(persistence)
import numpy.linalg as la
def logit_fit(X, y, it=50):
    beta = np.zeros(X.shape[1])
    for _ in range(it):
        eta = X @ beta; mu = 1 / (1 + np.exp(-eta)); W = mu * (1 - mu)
        H = X.T @ (X * W[:, None]); g = X.T @ (y - mu)
        step = la.solve(H, g); beta += step
        if np.max(np.abs(step)) < 1e-10: break
    return beta
def design(df):
    X = np.column_stack([np.ones(len(df))] + [(df.regime.values == k).astype(float) for k in (1, 2, 3)] + [np.log2(df.pers.values)])
    return X
beta = logit_fit(design(a), a.recurving.values.astype(float))[-1]
seasons = np.arange(1966, 2026); blocks = [seasons[i:i + 3] for i in range(0, 60, 3)]
early = [bl for bl in blocks if bl[0] <= 1995]; late = [bl for bl in blocks if bl[0] >= 1996]
bs = []
for _ in range(4999):
    pick = [early[i] for i in rng.integers(0, len(early), len(early))] + [late[i] for i in rng.integers(0, len(late), len(late))]
    idx = np.concatenate([np.flatnonzero(a.season.values == s) for bl in pick for s in bl])
    try: bs.append(logit_fit(design(a.iloc[idx]), a.recurving.values[idx].astype(float))[-1])
    except la.LinAlgError: pass
lo, hi = np.percentile(bs, [2.5, 97.5])
by = a.groupby(['regime', 'recurving']).pers.mean().unstack()
res = dict(n_storms=int(len(a)), persistence_days_median=float(np.median(p)), persistence_days_mean=float(p.mean()),
           left_censored_at_1_june=int(sum(1 for d, v in zip(a.adate, a.pers) if (pd.Timestamp(d) - pd.Timestamp(d.year, 6, 1)).days + 1 == v)),
           within_regime_mean_difference_days=obs, permutation_p=pval, n_permutations=B,
           mean_persistence_by_regime={f'regime_{k+1}': {'nonrecurving': float(by.loc[k, 0]), 'recurving': float(by.loc[k, 1])} for k in range(4)},
           logistic_or_per_doubling=float(np.exp(beta)), logistic_or_ci95=[float(np.exp(lo)), float(np.exp(hi))], n_bootstrap=len(bs))
json.dump(res, open(out, 'w'), indent=2); print(json.dumps(res, indent=1))
