'Exploratory annual-steering counterfactual analysis; not preregistered.\n\nEarly-period storm seeds are integrated with randomly sampled early-year monthly winds for group A and late-year monthly winds for group B. Both groups average N_REAL realizations under a common model calibration. The contrast B-A describes annual-steering distribution differences.\nMonthly fields do not resolve weather-scale steering changes; this analysis is not a six-hourly environmental-steering attribution experiment.'

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = _LAYOUT_ROOT / 'analysis/models'
EXP = _LAYOUT_ROOT

spec = importlib.util.spec_from_file_location("cfs", _LAYOUT_ROOT / 'analysis/models/counterfactual_steering.py')
cfs = importlib.util.module_from_spec(spec)
sys.modules["cfs"] = cfs
spec.loader.exec_module(cfs)

N_REAL = 20   # 随机年份配对的集合数
SEED = 20260804


def _interp2(field: np.ndarray, lat0: float, dlat: float, nlat: int,
             lon0: float, dlon: float, nlon: int, la: float, lo: float) -> float:
    'Scalar bilinear interpolation equivalent to cfs.bilinear.'
    fy = (la - lat0) / dlat
    fx = (lo - lon0) / dlon
    j = int(fy); i = int(fx)
    if j < 0 or j >= nlat - 1 or i < 0 or i >= nlon - 1:
        return float("nan")
    wy = fy - j; wx = fx - i
    return (field[j, i] * (1 - wy) * (1 - wx)
            + field[j, i + 1] * (1 - wy) * wx
            + field[j + 1, i] * wy * (1 - wx)
            + field[j + 1, i + 1] * wy * wx)


def integrate_yearly(seeds: pd.DataFrame, year_map: dict, cache: dict,
                     par: dict, region: dict) -> pd.DataFrame:
    """逐风暴使用其指定年份的月场做平流。"""
    rows = []
    a, ub, vb = par["alpha"], par["u_beta"], par["v_beta"]
    lat, lon = cache["lat"], cache["lon"]
    lat0, dlat, nlat = lat[0], lat[1] - lat[0], len(lat)
    lon0, dlon, nlon = lon[0], lon[1] - lon[0], len(lon)
    u_st, v_st = cache["u"], cache["v"]
    # (year, month) -> 月场行号，避免每个风暴重复线性搜索
    idx_of = {(y, m): k for k, (y, m) in enumerate(zip(cache["yrs"], cache["mons"]))}
    lo_min, lo_max = region["lon_min"], region["lon_max"]
    la_min, la_max = region["lat_min"], region["lat_max"]
    for r in seeds.itertuples():
        k = idx_of.get((year_map[r.sid], r.month))
        if k is None:
            continue
        u_f = u_st[k]; v_f = v_st[k]
        la, lo = r.lat, r.lon
        for _ in range(int(r.steps)):
            if not (lo_min <= lo <= lo_max and la_min <= la <= la_max):
                break
            rows.append((r.sid, r.season, lo, la))
            su = _interp2(u_f, lat0, dlat, nlat, lon0, dlon, nlon, la, lo)
            sv = _interp2(v_f, lat0, dlat, nlat, lon0, dlon, nlon, la, lo)
            if not (su == su and sv == sv):   # NaN 检查
                break
            lo += (a * su + ub) * cfs.DT / (cfs.EARTH_R * np.cos(la * cfs.RAD)) / cfs.RAD
            la += (a * sv + vb) * cfs.DT / cfs.EARTH_R / cfs.RAD
    return pd.DataFrame(rows, columns=["sid", "season", "lon", "lat"])


def main() -> None:
    cfg = cfs.load_config()
    tracks, region = cfs.load_tracks(cfg)
    lon_e, lat_e = cfs.grid_edges(cfg, region)
    early_years = np.arange(cfs.START, cfs.SPLIT + 1)
    late_years = np.arange(cfs.SPLIT + 1, cfs.END + 1)

    f1 = cfs.annual_relative_field(tracks[tracks.season <= cfs.SPLIT], early_years, lon_e, lat_e)
    f2 = cfs.annual_relative_field(tracks[tracks.season > cfs.SPLIT], late_years, lon_e, lat_e)
    d_obs = f2.mean(0) - f1.mean(0)

    cache = cfs._load_steering_once()
    # 直接读主试验存档的标定值：既省去重复标定的开销，也保证参数与主试验完全一致
    prev = json.loads((_LAYOUT_ROOT / 'results/models/counterfactual/criteria.json').read_text(encoding="utf-8"))
    par = prev["calibration"]
    print(f"沿用主试验标定 α={par['alpha']:.3f} beta=({par['u_beta']:+.2f},{par['v_beta']:+.2f})"
          f"（读自 results/criteria.json，未重新标定）")

    g = tracks.groupby("sid")
    seeds = g.first().reset_index()[["sid", "season", "lat", "lon", "iso_time"]]
    seeds["month"] = seeds["iso_time"].dt.month
    seeds["steps"] = g.size().values
    s1 = seeds[seeds.season <= cfs.SPLIT].reset_index(drop=True)

    # 对称设计：A、B 两组都用随机抽年并取 N_REAL 次平均。
    # 若 A 只用"各风暴自己那一年"（单次实现）而 B 取多次平均，两组平滑程度不同，
    # 差值会混入 A 的抽样噪声，使相关系数被压向 0。必须对称。
    rng = np.random.default_rng(SEED)

    def ensemble(years: np.ndarray, tag: str) -> tuple[np.ndarray, float]:
        acc = None
        reals = []
        for k in range(N_REAL):
            m = {r.sid: int(rng.choice(years)) for r in s1.itertuples()}
            syn = integrate_yearly(s1, m, cache, par, region)
            f = cfs.annual_relative_field(syn, early_years, lon_e, lat_e).mean(0)
            reals.append(f)
            acc = f if acc is None else acc + f
        stack = np.asarray(reals)
        # 实现间标准差的场均值：衡量该组自身的抽样噪声水平
        noise = float(stack.std(0).mean())
        print(f"  {tag} 完成 {N_REAL} 次实现，实现间噪声 {noise:.3e}")
        return acc / N_REAL, noise

    fld_a, noise_a = ensemble(early_years, "参照组A(随机早期年)")
    fld_b, noise_b = ensemble(late_years, "反事实B(随机晚期年)")

    d_ia = fld_b - fld_a
    r_ia = float(np.corrcoef(d_ia.ravel(), d_obs.ravel())[0, 1])
    ratio = cfs.tvd(d_ia) / cfs.tvd(d_obs)

    print("\n== 年际分布反事实（事后探索，非预注册）==")
    print(f"  r(d_ia, d_obs)          = {r_ia:+.3f}")
    print(f"  TVD(d_ia)/TVD(d_obs)    = {ratio:.3f}")
    print(f"  对照：气候态版 r        = {prev['r_steer_obs']:+.3f}，TVD比 {prev['tvd_ratio_steer']:.3f}")
    gain = "是" if abs(r_ia) > abs(prev["r_steer_obs"]) + 0.05 else "否"
    print(f"  换成年际分布后是否实质改善：{gain}")

    out = {
        "note": "事后探索，非预注册；不得作为预设判据结果引用",
        "n_realizations": N_REAL,
        "seed": SEED,
        "r_ia_obs": r_ia,
        "noise_a": noise_a,
        "noise_b": noise_b,
        "signal_to_noise": float(np.abs(d_ia).mean() / max(noise_a, noise_b)),
        "tvd_ratio_ia": ratio,
        "ref_r_steer_clim": prev["r_steer_obs"],
        "ref_tvd_ratio_clim": prev["tvd_ratio_steer"],
    }
    (_LAYOUT_ROOT / 'results/models/counterfactual/interannual.json').write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez(_LAYOUT_ROOT / 'results/models/counterfactual/interannual_fields.npz', d_ia=d_ia, d_obs=d_obs)
    print(f"\n结果写入 {_LAYOUT_ROOT / 'results/models/counterfactual/interannual.json'}")


if __name__ == "__main__":
    main()
