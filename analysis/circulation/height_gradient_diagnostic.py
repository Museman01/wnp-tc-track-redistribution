"""June-October 500-hPa height-gradient diagnostic for the subtropical high.

Tests whether the 1966-1995 to 1996-2025 height rise is spatially uniform or changes
the horizontal height gradients that set geostrophic steering.
Input: data/derived/circulation_monthly_coarse.npz (2.5 deg ERA5 monthly z500, m).
Usage: python -B height_gradient_diagnostic.py <npz> <out_json>
"""
import json, sys
import numpy as np

R = 6.371e6
SEED, NPERM, NBOOT, BLOCK = 202406, 9999, 4999, 3


def season_mean(d):
    yrs = np.unique(d["year"]); z = d["z500"]; m = d["month"]
    out = np.stack([z[(d["year"] == y) & (m >= 6) & (m <= 10)].mean(0) for y in yrs])
    return yrs, out


def gradients(Z, lat, lon):
    """Centered differences in m per 1000 km. lat is descending."""
    phi = np.deg2rad(lat)[:, None]
    dy = np.gradient(np.deg2rad(lat)) * R
    dx = np.gradient(np.deg2rad(lon)) * R
    gy = np.gradient(Z, axis=-2) / dy[:, None] * 1e6          # dZ/dy (northward positive)
    gx = np.gradient(Z, axis=-1) / (dx[None, :] * np.cos(phi)) * 1e6
    return gx, gy


def box_mean(F, lat, lon, la, lo):
    sl = (lat >= la[0]) & (lat <= la[1]); so = (lon >= lo[0]) & (lon <= lo[1])
    w = np.cos(np.deg2rad(lat[sl]))[:, None] * np.ones(so.sum())[None, :]
    sub = F[..., sl, :][..., so]
    return (sub * w).sum((-2, -1)) / w.sum()


def block_index(n, rng):
    starts = rng.integers(0, n - BLOCK + 1, size=int(np.ceil(n / BLOCK)))
    return np.concatenate([np.arange(s, s + BLOCK) for s in starts])[:n]


def period_stats(series, early, late, rng):
    a, b = series[early], series[late]
    diff = b.mean() - a.mean()
    boots = np.array([b[block_index(len(b), rng)].mean() - a[block_index(len(a), rng)].mean()
                      for _ in range(NBOOT)])
    allv = np.concatenate([a, b]); n = len(allv)
    blocks = [allv[i:i + BLOCK] for i in range(0, n, BLOCK)]
    na = len(a) // BLOCK; cnt = 0
    for _ in range(NPERM):
        order = rng.permutation(len(blocks))
        x = np.concatenate([blocks[i] for i in order])
        cnt += abs(x[na * BLOCK:].mean() - x[:na * BLOCK].mean()) >= abs(diff) - 1e-12
    return {"early_mean": float(a.mean()), "late_mean": float(b.mean()), "difference": float(diff),
            "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            "perm_p": float((cnt + 1) / (NPERM + 1))}


def main(npz, out):
    d = np.load(npz)
    lat, lon = d["latitude"], d["longitude"]
    yrs, Z = season_mean(d)
    early = (yrs >= 1966) & (yrs <= 1995); late = (yrs >= 1996) & (yrs <= 2025)
    gx, gy = gradients(Z, lat, lon)
    rng = np.random.default_rng(SEED)
    boxes = {
        "domain_110_180E_10_40N": ((10, 40), (110, 180)),
        "southern_flank_10_20N_120_160E": ((10, 20), (120, 160)),
        "western_flank_15_30N_110_130E": ((15, 30), (110, 130)),
        "northern_side_25_40N_110_180E": ((25, 40), (110, 180)),
    }
    res = {"input": npz, "season": "June-October", "years": [int(yrs[0]), int(yrs[-1])],
           "grid_deg": float(abs(lat[1] - lat[0])), "units": "Z in m; gradients in m per 1000 km",
           "seed": SEED, "n_perm": NPERM, "n_boot": NBOOT, "block_years": BLOCK}
    dZ = Z[late].mean(0) - Z[early].mean(0)
    la, lo = boxes["domain_110_180E_10_40N"]
    sl = (lat >= la[0]) & (lat <= la[1]); so = (lon >= lo[0]) & (lon <= lo[1])
    w = np.cos(np.deg2rad(lat[sl]))[:, None] * np.ones(so.sum())[None, :]
    sub = dZ[sl][:, so]; mu = float((sub * w).sum() / w.sum())
    sd = float(np.sqrt(((sub - mu) ** 2 * w).sum() / w.sum()))
    res["height_change_field"] = {"area_mean_m": mu, "spatial_sd_m": sd, "min_m": float(sub.min()), "max_m": float(sub.max())}
    for name, (la, lo) in boxes.items():
        r = {"Z": period_stats(box_mean(Z, lat, lon, la, lo), early, late, rng)}
        if name.startswith(("southern", "northern")):
            r["dZdy"] = period_stats(box_mean(gy, lat, lon, la, lo), early, late, rng)
        if name.startswith("western"):
            r["dZdx"] = period_stats(box_mean(gx, lat, lon, la, lo), early, late, rng)
        res[name] = r
    json.dump(res, open(out, "w"), indent=2)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
