'Storm time sampling.'


from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT



import json

import sys




import netCDF4

import numpy as np

import pandas as pd

import yaml



HERE = _LAYOUT_ROOT / 'analysis'

PROJECT = _LAYOUT_ROOT

TASK = PROJECT

sys.path.insert(0, str(PROJECT))



from analysis.shared.agency_data import read_ibtracs_agencies, build_agency_catalog  # noqa: E402



AGENCIES = ["USA", "TOKYO", "CMA"]

AGENCY_LABEL = {"USA": "USA", "TOKYO": "JMA", "CMA": "CMA"}

P1 = (1966, 1995)

P2 = (1996, 2025)





def load_config() -> dict:

    with open(PROJECT / "config" / "config.yaml", "r", encoding="utf-8") as fh:

        return yaml.safe_load(fh)





# --------------------------------------------------------------------------

# 场：把 0.25° 月场按 cos(lat) 权重聚合到 2.5° 分析格

# --------------------------------------------------------------------------

def coarsen_steering(cfg: dict, months, years) -> dict:

    region = cfg["regions"]["tc"]

    bw = cfg["grids"]["density_bin"]

    lon_edges = np.arange(region["lon_min"], region["lon_max"] + bw, bw)

    lat_edges = np.arange(region["lat_min"], region["lat_max"] + bw, bw)



    # netCDF4 直读：时间维按块连续读取。xarray 的布尔时间选择会对压缩块做

    # 非连续花式索引，导致同一块被反复解压，这里避免该路径。

    nc = netCDF4.Dataset(PROJECT / "data" / "interim" / "steering.nc")

    all_times = pd.to_datetime(

        netCDF4.num2date(

            nc.variables["time"][:],

            nc.variables["time"].units,

            nc.variables["time"].calendar,

            only_use_cftime_datetimes=False,

        )

    )

    lat_full = nc.variables["latitude"][:].astype(float)

    lon_full = nc.variables["longitude"][:].astype(float)

    # 文件内纬度递减；取窗口后统一转为递增

    lat_sel = np.where(

        (lat_full >= region["lat_min"]) & (lat_full <= region["lat_max"])

    )[0]

    lon_sel = np.where(

        (lon_full >= region["lon_min"]) & (lon_full <= region["lon_max"])

    )[0]

    la0, la1 = lat_sel.min(), lat_sel.max() + 1

    lo0, lo1 = lon_sel.min(), lon_sel.max() + 1

    lat = lat_full[la0:la1]

    lon = lon_full[lo0:lo1]

    flip = lat[0] > lat[-1]

    if flip:

        lat = lat[::-1]



    keep = np.where(

        np.isin(np.asarray(all_times.month), months)

        & np.isin(np.asarray(all_times.year), years)

    )[0]

    times = all_times[keep]

    w = np.cos(np.deg2rad(lat))



    nlat, nlon = len(lat_edges) - 1, len(lon_edges) - 1

    jlat = np.clip(np.digitize(lat, lat_edges, right=False) - 1, 0, nlat - 1)

    ilon = np.clip(np.digitize(lon, lon_edges, right=False) - 1, 0, nlon - 1)



    # 聚合矩阵：A_lat (nlat, nlat_fine) 带 cos(lat) 权重；A_lon (nlon, nlon_fine) 等权

    A_lat = np.zeros((nlat, len(lat)))

    A_lat[jlat, np.arange(len(lat))] = w

    A_lon = np.zeros((nlon, len(lon)))

    A_lon[ilon, np.arange(len(lon))] = 1.0

    den = A_lat.sum(axis=1)[:, None] * A_lon.sum(axis=1)[None, :]



    def agg(name: str) -> np.ndarray:

        var = nc.variables[name]

        tchunk = var.chunking()[0] if isinstance(var.chunking(), list) else 344

        out = np.empty((len(keep), nlat, nlon))

        for s in range(0, var.shape[0], tchunk):

            e = min(s + tchunk, var.shape[0])

            sel = keep[(keep >= s) & (keep < e)]

            if sel.size == 0:

                continue

            block = np.asarray(var[s:e, la0:la1, lo0:lo1], dtype=float)

            if flip:

                block = block[:, ::-1, :]

            sub = block[sel - s]

            num = np.einsum("jJ,tJI,iI->tji", A_lat, sub, A_lon, optimize=True)

            pos = np.searchsorted(keep, sel)

            out[pos] = num / np.where(den == 0, np.nan, den)

            del block, sub

        return out



    out = {

        "u": agg("u_steer"),

        "v": agg("v_steer"),

        "year": times.year.values,

        "month": times.month.values,

        "lon_edges": lon_edges,

        "lat_edges": lat_edges,

    }

    nc.close()

    return out





def coarsen_steering_cached(cfg, months, years):

    '将月引导风聚合到分析网格，并保存共享处理场。'

    cache = _LAYOUT_ROOT / 'results/circulation/storm_sampling' / "steering_2p5deg_cache.npz"

    cache.parent.mkdir(parents=True, exist_ok=True)

    if cache.exists():

        z = np.load(cache)

        return {k: z[k] for k in z.files}

    grid = coarsen_steering(cfg, months, years)

    np.savez_compressed(cache, **grid)

    return grid





# --------------------------------------------------------------------------

# 分布：逐年归一的路径点分布 n_y(month, lat_bin, lon_bin)

# --------------------------------------------------------------------------

def annual_point_distributions(tracks: pd.DataFrame, cfg: dict, months, years, grid):

    region = cfg["regions"]["tc"]

    lon_edges, lat_edges = grid["lon_edges"], grid["lat_edges"]

    nm, nlat, nlon = len(months), len(lat_edges) - 1, len(lon_edges) - 1



    t = tracks.copy()

    t["iso_time"] = pd.to_datetime(t["iso_time"])

    t["month"] = t["iso_time"].dt.month

    t = t.loc[

        t["month"].isin(months)

        & t["season"].isin(years)

        & t["lon"].between(region["lon_min"], region["lon_max"])

        & t["lat"].between(region["lat_min"], region["lat_max"])

    ]



    month_index = {m: k for k, m in enumerate(months)}

    dists, kept_years, npts = {}, [], {}

    for y in years:

        g = t.loc[t["season"] == y]

        if len(g) == 0:

            continue

        arr = np.zeros((nm, nlat, nlon))

        mk = g["month"].map(month_index).to_numpy()

        jj = np.clip(np.digitize(g["lat"].to_numpy(), lat_edges) - 1, 0, nlat - 1)

        ii = np.clip(np.digitize(g["lon"].to_numpy(), lon_edges) - 1, 0, nlon - 1)

        np.add.at(arr, (mk, jj, ii), 1.0)

        total = arr.sum()

        if total == 0:

            continue

        dists[y] = arr / total          # 年度归一化：每年权重相等

        npts[y] = int(total)

        kept_years.append(y)

    return dists, kept_years, npts





def field_by_year_month(grid, months, years):

    """U[year_index, month_index, lat, lon]，缺格保留 NaN。"""

    nm = len(months)

    nlat = len(grid["lat_edges"]) - 1

    nlon = len(grid["lon_edges"]) - 1

    yindex = {y: k for k, y in enumerate(years)}

    mindex = {m: k for k, m in enumerate(months)}

    U = {c: np.full((len(years), nm, nlat, nlon), np.nan) for c in ("u", "v")}

    for k, (y, m) in enumerate(zip(grid["year"], grid["month"])):

        if y not in yindex or m not in mindex:

            continue

        for c in ("u", "v"):

            U[c][yindex[y], mindex[m]] = grid[c][k]

    return U





# --------------------------------------------------------------------------

# 分解

# --------------------------------------------------------------------------

def decompose(U, dists, years_all, years_p1, years_p2, valid):

    """返回一个分量（u 或 v）的三项分解与固定空间版本。"""

    yindex = {y: k for k, y in enumerate(years_all)}



    def stack(years):

        n = np.array([dists[y] for y in years])            # (ny, m, j, i)

        f = np.array([U[yindex[y]] for y in years])        # (ny, m, j, i)

        return n, f



    n1, U1 = stack(years_p1)

    n2, U2 = stack(years_p2)



    # 有效格（两期都非 NaN），并把权重重归一到有效格上

    def renorm(n):

        m = n * valid

        s = m.sum(axis=(1, 2, 3), keepdims=True)

        return m / np.where(s == 0, np.nan, s)



    n1, n2 = renorm(n1), renorm(n2)

    U1 = np.where(valid, U1, 0.0)

    U2 = np.where(valid, U2, 0.0)



    f1, f2 = n1.mean(axis=0), n2.mean(axis=0)

    Uc1, Uc2 = U1.mean(axis=0), U2.mean(axis=0)



    S1_clim = float((f1 * Uc1).sum())

    S2_clim = float((f2 * Uc2).sum())

    S1_act = float(np.mean([(n1[k] * U1[k]).sum() for k in range(len(years_p1))]))

    S2_act = float(np.mean([(n2[k] * U2[k]).sum() for k in range(len(years_p2))]))



    dU = Uc2 - Uc1

    df = f2 - f1

    field = float((0.5 * (f1 + f2) * dU).sum())

    sampling = float((0.5 * df * (Uc1 + Uc2)).sum())

    cosampling = (S2_act - S2_clim) - (S1_act - S1_clim)



    # ---- 空间固定版本：空间边缘锁死为两期合并分布，只让月组成与年份变化 ----

    def split_marginals(f):

        month_share = f.sum(axis=(1, 2))                  # (m,)

        with np.errstate(invalid="ignore", divide="ignore"):

            spatial = f / month_share[:, None, None]

        return month_share, np.nan_to_num(spatial)



    ms1, _ = split_marginals(f1)

    ms2, _ = split_marginals(f2)

    f_pool = 0.5 * (f1 + f2)

    _, sp_pool = split_marginals(f_pool)

    g1 = ms1[:, None, None] * sp_pool                     # 固定空间、用 P1 月组成

    g2 = ms2[:, None, None] * sp_pool                     # 固定空间、用 P2 月组成

    seasonal = float((0.5 * (g2 - g1) * (Uc1 + Uc2)).sum())

    field_fixed = float((0.5 * (g1 + g2) * dU).sum())



    # 固定空间的年际共采样：用 g_p 权重取实际年月场

    def cosamp_fixed(n_stack, U_stack, g):

        ms = np.array([nk.sum(axis=(1, 2)) for nk in n_stack])   # (ny, m) 各年月组成

        acts = []

        for k in range(len(n_stack)):

            gk = ms[k][:, None, None] * sp_pool

            acts.append((gk * U_stack[k]).sum())

        return float(np.mean(acts) - (g * np.mean(U_stack, axis=0)).sum())



    cos1 = cosamp_fixed(n1, U1, g1)

    cos2 = cosamp_fixed(n2, U2, g2)



    return {

        "S1_actual": S1_act,

        "S2_actual": S2_act,

        "S1_clim": S1_clim,

        "S2_clim": S2_clim,

        "total_change": S2_act - S1_act,

        "field": field,

        "sampling": sampling,

        "cosampling": cosampling,

        "additivity_residual": (S2_act - S1_act) - (field + sampling + cosampling),

        "fixed_space": {

            "field": field_fixed,

            "seasonal": seasonal,

            "cosampling": cos2 - cos1,

        },

    }





def block_bootstrap_years(years, block, rng):

    years = np.asarray(years)

    n = len(years)

    starts = rng.integers(0, n, size=int(np.ceil(n / block)))

    out = np.concatenate([years[np.arange(s, s + block) % n] for s in starts])[:n]

    return out





def main():

    cfg = load_config()

    months = cfg["typhoon_season"]

    years_all = list(range(P1[0], P2[1] + 1))

    years_p1 = [y for y in years_all if P1[0] <= y <= P1[1]]

    years_p2 = [y for y in years_all if P2[0] <= y <= P2[1]]



    grid = coarsen_steering_cached(cfg, months, years_all)

    U = field_by_year_month(grid, months, years_all)

    valid = np.isfinite(U["u"]).all(axis=0) & np.isfinite(U["v"]).all(axis=0)

    Uf = {c: np.nan_to_num(U[c]) for c in ("u", "v")}



    # σ 参考：主季月际标准差（有效格、按 cos 权重的空间平均后取时间标准差）

    frame = read_ibtracs_agencies(

        PROJECT / "data" / "raw" / "ibtracs_wp.csv",

        start=P1[0], end=P2[1],

    )



    out = {

        "note": "事后探索，非预注册；不得作为预设判据结果引用",

        "periods": {"P1": list(P1), "P2": list(P2)},

        "months": months,

        "bin_deg": cfg["grids"]["density_bin"],

        "region": cfg["regions"]["tc"],

        "n_valid_cells": int(valid.sum()),

        "bootstrap": {

            "samples": cfg["statistics"]["bootstrap_samples"],

            "block": cfg["statistics"]["bootstrap_block"],

            "seed": cfg["statistics"]["random_seed"],

        },

        "agencies": {},

    }



    nboot = cfg["statistics"]["bootstrap_samples"]

    block = cfg["statistics"]["bootstrap_block"]



    

    # 两期平均场差与月际标准差之比。

    yidx = {y: k for k, y in enumerate(years_all)}

    i1 = [yidx[y] for y in years_p1]

    i2 = [yidx[y] for y in years_p2]

    ref = {}

    for comp in ("u", "v"):

        arr = Uf[comp]                                    # (ny, m, j, i)

        area = arr.mean(axis=(2, 3))                      # 域平均，逐年逐月

        d = float(area[i2].mean() - area[i1].mean())

        sd = float(area.reshape(-1).std(ddof=1))

        ref[comp] = {

            "period_mean_difference_ms": d,

            "month_to_month_sd_ms": sd,

            "delta_over_sigma": abs(d) / sd,

        }

    out["reference_diagnostic_area_mean"] = ref



    for agency in AGENCIES:

        cat = build_agency_catalog(frame, agency)

        tracks = cat["tracks"]

        dists, kept, npts = annual_point_distributions(

            tracks, cfg, months, years_all, grid

        )

        yp1 = [y for y in years_p1 if y in dists]

        yp2 = [y for y in years_p2 if y in dists]



        res = {}

        for comp in ("u", "v"):

            res[comp] = decompose(Uf[comp], dists, years_all, yp1, yp2, valid)



        # 自助：期内 3 年块重采样

        rng = np.random.default_rng(cfg["statistics"]["random_seed"])

        boot = {c: {k: [] for k in ("total_change", "field", "sampling", "cosampling",

                                    "seasonal_fixed", "cosampling_fixed")}

                for c in ("u", "v")}

        for _ in range(nboot):

            b1 = list(block_bootstrap_years(yp1, block, rng))

            b2 = list(block_bootstrap_years(yp2, block, rng))

            for comp in ("u", "v"):

                d = decompose(Uf[comp], dists, years_all, b1, b2, valid)

                boot[comp]["total_change"].append(d["total_change"])

                boot[comp]["field"].append(d["field"])

                boot[comp]["sampling"].append(d["sampling"])

                boot[comp]["cosampling"].append(d["cosampling"])

                boot[comp]["seasonal_fixed"].append(d["fixed_space"]["seasonal"])

                boot[comp]["cosampling_fixed"].append(d["fixed_space"]["cosampling"])



        ci = {}

        for comp in ("u", "v"):

            ci[comp] = {}

            for k, vals in boot[comp].items():

                a = np.asarray(vals)

                ci[comp][k] = {

                    "p2.5": float(np.percentile(a, 2.5)),

                    "p97.5": float(np.percentile(a, 97.5)),

                    "sd": float(a.std(ddof=1)),

                    "frac_same_sign": float(np.mean(np.sign(a) == np.sign(a.mean()))),

                }



        out["agencies"][AGENCY_LABEL[agency]] = {

            "n_years_P1": len(yp1),

            "n_years_P2": len(yp2),

            "n_points_P1": int(sum(npts[y] for y in yp1)),

            "n_points_P2": int(sum(npts[y] for y in yp2)),

            "decomposition": res,

            "bootstrap_ci": ci,

        }



    dest = _LAYOUT_ROOT / 'results/circulation/storm_sampling'

    dest.mkdir(parents=True, exist_ok=True)

    with open(dest / "storm_time_sampling.json", "w", encoding="utf-8") as fh:

        json.dump(out, fh, indent=2, ensure_ascii=False)

    print(json.dumps(out["agencies"], indent=2, ensure_ascii=False)[:4000])





if __name__ == "__main__":

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    main()

