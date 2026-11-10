"""Storm-matched wind statistics and circulation-conditioned path decomposition."""
import numpy as np
import pandas as pd
from analysis.shared.block_resampling import block_order, moving_block_sample, bh_fdr
from analysis.tracks.pathway_closure import decompose_period_means, block_indices
from analysis.shared.inference_stats import compositional_change_test
from analysis.shared.agency_catalogues import annual_fields
YEARS=np.arange(1966,2026)
EARLY=np.arange(30)
LATE=np.arange(30,60)
SEED=202406
GEOMETRIES={'annulus_5_7':(5.,7.),'disk_5':(0.,5.),'annulus_3_7':(3.,7.)}
METRICS=['corr_u','corr_v','corr_vector','residual_u','residual_v','rmse_vector']


def motion_pairs(df):
    out = []
    for sid, g in df.sort_values(['sid', 'iso_time']).groupby('sid'):
        for a, b in zip(g.itertuples(), list(g.itertuples())[1:]):
            dt = (b.iso_time - a.iso_time).total_seconds()
            if 0 < dt <= 43200 and 100 <= a.lon <= 180 and (0 <= a.lat <= 40):
                dl = (b.lon - a.lon + 180) % 360 - 180
                out.append(dict(sid=sid, season=a.season, iso_time=a.iso_time, lat=a.lat, lon=a.lon, obs_u=6371000 * np.deg2rad(dl) * np.cos(np.deg2rad((a.lat + b.lat) / 2)) / dt, obs_v=6371000 * np.deg2rad(b.lat - a.lat) / dt, dt_hours=dt / 3600))
    return pd.DataFrame(out)

def moment(a, b):
    return np.array([len(a), a.sum(), b.sum(), (a * a).sum(), (b * b).sum(), (a * b).sum()], float)

def corr(m):
    n, sx, sy, xx, yy, xy = m
    den = (xx - sx * sx / n) * (yy - sy * sy / n)
    return (xy - sx * sy / n) / np.sqrt(den) if n > 1 and den > 0 else np.nan

def vector_corr(m):
    n = m[:, 0]
    cov = m[:, 5] - m[:, 1] * m[:, 2] / n
    va = m[:, 3] - m[:, 1] ** 2 / n
    vb = m[:, 4] - m[:, 2] ** 2 / n
    den = va.sum() * vb.sum()
    return cov.sum() / np.sqrt(den) if np.all(n > 1) and den > 0 else np.nan

def regime_annual(frame, k):
    expected = pd.DatetimeIndex(np.concatenate([pd.date_range(f'{y}-06-01', f'{y}-10-31').values for y in YEARS]))
    if not np.array_equal(pd.to_datetime(frame.date).values, expected.values):
        raise ValueError('Daily regime axis must be complete sorted Jun-Oct, all years')
    rows = []
    for y, g in frame.groupby('year', sort=True):
        labels = g.regime.to_numpy()
        row = dict(year=int(y), ridge_latitude=g.ridge_latitude.mean(), favorable_fraction=g.favorable.mean())
        starts = np.r_[0, np.flatnonzero(labels[1:] != labels[:-1]) + 1]
        ends = np.r_[starts[1:], len(labels)]
        for c in range(k):
            ix = np.flatnonzero(labels[starts] == c)
            row[f'frequency_{c}'] = np.mean(labels == c)
            row[f'duration_{c}'] = float(np.mean((ends - starts)[ix])) if len(ix) else np.nan
            row[f'censored_runs_{c}'] = int(np.sum((starts[ix] == 0) | (ends[ix] == len(labels))))
        rows.append(row)
    return pd.DataFrame(rows)

def geometry_weights(lat, lon, center_lat, center_lon, geometry):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if geometry not in GEOMETRIES or lat.ndim != 1 or lon.ndim != 1 or (len(lat) < 2) or (len(lon) < 2):
        raise ValueError('Invalid geometry/grid')
    if not np.allclose(np.diff(lat), -1) or not np.allclose(np.diff(lon), 1):
        raise ValueError('Expected descending 1-degree latitude / ascending longitude')
    inner, outer = GEOMETRIES[geometry]
    if not np.isfinite([center_lat, center_lon]).all() or abs(center_lat) + outer >= 90:
        raise ValueError('Invalid center or polar ring')
    half_lon = np.rad2deg(np.arcsin(np.sin(np.deg2rad(outer)) / np.cos(np.deg2rad(center_lat))))
    covered = center_lat - outer >= lat.min() - 0.5 and center_lat + outer <= lat.max() + 0.5 and (center_lon - half_lon >= lon.min() - 0.5) and (center_lon + half_lon <= lon.max() + 0.5)
    if not covered:
        return (None, {'complete_geometry': False, 'geometric_longitude_half_width_deg': float(half_lon), 'selected_grid_cells': 0})
    a, b = (np.deg2rad(lat[:, None]), np.deg2rad(lon[None, :]))
    c, d = np.deg2rad([center_lat, center_lon])
    distance = np.rad2deg(np.arccos(np.clip(np.sin(a) * np.sin(c) + np.cos(a) * np.cos(c) * np.cos(b - d), -1, 1)))
    mask = (distance >= inner) & (distance <= outer)
    weights = mask * np.cos(a)
    if weights.sum() <= 0:
        raise ValueError('Empty geometry')
    return (weights / weights.sum(), {'complete_geometry': True, 'geometric_longitude_half_width_deg': float(half_lon), 'selected_grid_cells': int(mask.sum())})

def sampled_month_mean(ds, year, month):
    start = pd.Timestamp(year=year, month=month, day=1)
    end = start + pd.offsets.MonthBegin(1)
    times = pd.date_range(start, end, freq='6h', inclusive='left')
    fields = ['u_deep', 'v_deep', 'u_shallow', 'v_shallow']
    part = ds[fields].sel(time=times)
    if not np.array_equal(part.time.values, times.values):
        raise ValueError('Incomplete month timestamp axis')
    sums = {f: np.zeros((ds.sizes['latitude'], ds.sizes['longitude']), dtype='f8') for f in fields}
    for start in range(0, len(times), 16):
        for f in fields:
            v = part[f].isel(time=slice(start, start + 16)).values
            if not np.isfinite(v).all():
                raise ValueError('Nonfinite monthly source')
            sums[f] += v.sum(axis=0, dtype='f8')
    return {f: a / len(times) for f, a in sums.items()}

def extract_point(ds, instant, lat, lon, monthly_fields):
    point = ds.sel(time=pd.Timestamp(instant))
    out = {}
    for geometry in GEOMETRIES:
        weight, qc = geometry_weights(ds.latitude.values, ds.longitude.values, lat, lon, geometry)
        out[geometry + '_coverage'] = qc['complete_geometry']
        out[geometry + '_grid_cells'] = qc['selected_grid_cells']
        for layer in ('deep', 'shallow'):
            for component in 'uv':
                field = component + '_' + layer
                name = f'{geometry}_{layer}_{component}'
                if weight is None:
                    out[name] = np.nan
                    out[name + '_monthmatched'] = np.nan
                    continue
                a = point[field].values
                if not np.isfinite(a).all():
                    raise ValueError('Nonfinite instantaneous field')
                out[name] = float(np.sum(a * weight))
                out[name + '_monthmatched'] = float(np.sum(monthly_fields[field] * weight))
    return out

def stratified_44_order(rng):
    triples = [np.arange(i, i + 3) for start in (0, 22) for i in range(start, start + 21, 3)]
    order = rng.permutation(14)
    tails = rng.permutation([21, 43])
    return np.r_[np.concatenate([triples[i] for i in order[:7]]), tails[0], np.concatenate([triples[i] for i in order[7:]]), tails[1]]

def annual_sufficient(frame, columns, years):
    moments = []
    residuals = []
    for year in years:
        g = frame[frame.season == year]
        a = g[columns].to_numpy(float)
        b = g[['obs_u', 'obs_v']].to_numpy(float)
        if not len(g) or not np.isfinite(np.c_[a, b]).all():
            raise ValueError('Empty/nonfinite storm-year movement sample')
        moments.append([moment(a[:, j], b[:, j]) for j in range(2)])
        residuals.append([len(g), *(b - a).sum(axis=0), np.square(b - a).sum()])
    return (np.array(moments), np.array(residuals))

def pooled_metric(moments, residuals, indices):
    m = moments[indices].sum(axis=0)
    r = residuals[indices].sum(axis=0)
    return np.array([corr(m[0]), corr(m[1]), vector_corr(m), r[1] / r[0], r[2] / r[0], np.sqrt(r[3] / r[0])])

def infer(metric, years, nperm=9999, nboot=4999, guard=None, undefined_bootstrap='stop'):
    if undefined_bootstrap not in ('stop', 'report_na'):
        raise ValueError('Unknown undefined-CI policy')
    years = np.asarray(years)
    if not (np.array_equal(years, np.arange(1966, 2026)) or np.array_equal(years, np.arange(1982, 2026))):
        raise ValueError('Only locked60/44-year windows')
    n = len(years)
    half = n // 2
    early = np.arange(half)
    late = np.arange(half, n)
    e = metric(early)
    l = metric(late)
    delta = l - e
    if not np.isfinite(np.r_[e, l, delta]).all():
        raise ValueError('Undefined observed metric')
    exceed = np.zeros(len(delta), dtype=int)
    rng = np.random.default_rng(SEED)
    for iteration in range(nperm):
        if guard is not None and iteration % 128 == 0:
            guard()
        order = block_order(n, 3, rng) if n == 60 else stratified_44_order(rng)
        null = metric(order[half:]) - metric(order[:half])
        if not np.isfinite(null).all():
            raise ValueError('Undefined permutation metric')
        exceed += np.abs(null) >= np.abs(delta)
    rng = np.random.default_rng(SEED)
    samples = []
    for iteration in range(nboot):
        if guard is not None and iteration % 128 == 0:
            guard()
        a = moving_block_sample(early, 3, rng).astype(int)
        b = moving_block_sample(late, 3, rng).astype(int)
        value = metric(b) - metric(a)
        if not np.isfinite(value).all() and undefined_bootstrap == 'stop':
            raise ValueError('Undefined bootstrap metric')
        samples.append(value)
    samples = np.array(samples)
    bad = (~np.isfinite(samples)).sum(axis=0)
    ci = np.full((2, len(delta)), np.nan)
    for j in range(len(delta)):
        if bad[j] == 0:
            ci[:, j] = np.percentile(samples[:, j], [2.5, 97.5])
    return {'early': e, 'late': l, 'difference': delta, 'p': (exceed + 1) / (nperm + 1), 'ci_low': ci[0], 'ci_high': ci[1], 'nperm': nperm, 'nboot': nboot, 'window_years': n, 'bootstrap_undefined_counts': bad, 'bootstrap_policy': undefined_bootstrap}



def decompose(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.ndim == 2:
        a = a[:, :, None]
        b = b[:, :, None]
    p = a.sum((1, 2))
    q = b.sum((1, 2))
    if np.any(p <= 0) or np.any(q <= 0):
        raise ValueError('undefined conditional: zero regime support')
    delta = b.sum(0) - a.sum(0)
    if np.vdot(delta, delta) > 0:
        return decompose_period_means(p, q, a, b)
    f = a / p[:, None, None]
    g = b / q[:, None, None]
    between = ((q - p)[:, None, None] * (f + g) / 2).sum(0)
    within = ((q + p)[:, None, None] * (g - f) / 2).sum(0)
    return dict(delta=delta, between=between, within=within, pi_early=p, pi_late=q, f_early=f, f_late=g, closure_max_abs=float(abs(delta - between - within).max()))

def project(x, reference):
    den = float(np.vdot(reference, reference))
    return float(100 * np.vdot(x, reference) / den) if den > 1e-24 else np.nan

def permutation_order(n, rng):
    if n % 3:
        raise ValueError('year count must be divisible by3')
    return block_order(n, 3, rng)

def inference(j, nboot=4999, nperm=9999, seed=SEED):
    if j.shape[0] != 60 or not np.allclose(j.reshape(60, -1).sum(1), 1):
        raise ValueError('require60 annual unit joint fields')
    obs = decompose(j[EARLY].mean(0), j[LATE].mean(0))
    ref = obs['delta']
    rng = np.random.default_rng(seed)
    boot = np.empty((nboot, 2))
    for i in range(nboot):
        a = block_indices(EARLY, 3, rng)
        b = block_indices(LATE, 3, rng)
        r = decompose(j[a].mean(0), j[b].mean(0))
        boot[i] = [project(r[k], ref) for k in ['between', 'within']]
    rng = np.random.default_rng(seed)
    exceed = np.zeros(3, int)
    keys = ['between', 'within', 'delta']
    norms = np.array([np.linalg.norm(obs[k]) for k in keys])
    for i in range(nperm):
        o = permutation_order(60, rng)
        r = decompose(j[o[EARLY]].mean(0), j[o[LATE]].mean(0))
        exceed += np.array([np.linalg.norm(r[k]) for k in keys]) >= norms
    tv = compositional_change_test(j.sum(1).reshape(60, -1), EARLY, LATE, nperm=nperm, block=3, seed=seed)
    rows = []
    for i, k in enumerate(keys[:2]):
        ci = np.quantile(boot[:, i], [0.025, 0.975]) if np.isfinite(boot[:, i]).all() else [np.nan, np.nan]
        rows.append(dict(component=k, projection_percent=project(obs[k], ref), ci_low=ci[0], ci_high=ci[1], permutation_p=(exceed[i] + 1) / (nperm + 1), permutation_statistic='L2_norm_of_component_not_selected_projection'))
    return (obs, pd.DataFrame(rows), boot, tv, (exceed[2] + 1) / (nperm + 1))

def label_points(tracks, labels):
    t = tracks.copy()
    t['iso_time'] = pd.to_datetime(t.iso_time, utc=True)
    t = t[t.iso_time.dt.month.between(6, 10) & t.lon.between(100, 180) & t.lat.between(0, 40)].copy()
    t['regime'] = t.iso_time.dt.normalize().map(labels)
    if t.regime.isna().any():
        raise ValueError('missing point regime')
    return t

def path_joint(t, keys):
    j = np.zeros((60, len(keys), 16, 32))
    counts = []
    for i, y in enumerate(YEARS):
        g = t[t.season == y]
        n = len(g)
        if not n:
            raise ValueError(f'empty path year{y}')
        for k, r in enumerate(keys):
            h = g[g.regime == r]
            j[i, k] = np.histogram2d(h.lat, h.lon, bins=[np.arange(0, 40.1, 2.5), np.arange(100, 180.1, 2.5)])[0] / n
        counts.append(dict(year=y, points=n, storms=g.sid.nunique()))
    np.testing.assert_allclose(j.sum(1).reshape(60, -1), annual_fields(t, YEARS), atol=1e-14)
    return (j, pd.DataFrame(counts))
