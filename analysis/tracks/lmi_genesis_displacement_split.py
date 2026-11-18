"""Song & Klotzbach (2018)-style split of LMI latitude into genesis latitude and the
meridional displacement between genesis and LMI, applied to the study records.

Input : native tropical-storm-qualified six-hourly reports and native winds from
        results/ocean/ohc300/tracks_agency_definitions.csv.gz for PRIMARY, USA,
        JMA, and CMA. HKO uses the native tropical-storm points from the agency catalogue.
Rules : domain 100-180E, 0-40N (same as the track-composition analysis);
        genesis = first qualified position in the domain (operational genesis, Section 2.4);
        LMI = position of the maximum native wind in the domain; tied maxima are resolved
        three ways (first, last, mean latitude of the tied positions);
        storms are averaged within each year, then years are averaged (annual-mean basis,
        as in Song & Klotzbach 2018).
Output: results/tracks/lmi_genesis_displacement_split.csv
        two-period difference (1996-2025 minus 1966-1995) and least-squares linear trends of
        annual means over 1966-2016 and 1966-2025; displacement share = d(disp)/d(LMI lat).
Usage : python -m analysis.tracks.lmi_genesis_displacement_split <track CSV> <new output CSV>
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.shared.agency_catalogues import all_catalogues

def storms(path):
    tracks = pd.read_csv(path, parse_dates=["iso_time"])
    native = tracks.loc[tracks.definition.eq("ts_only")].copy()
    if set(native.agency) != {"PRIMARY", "USA", "JMA", "CMA"}:
        raise ValueError("The track CSV must contain PRIMARY, USA, JMA, and CMA.")
    hko = all_catalogues(1966)["HKO"]["ts_points"].assign(agency="HKO", definition="ts_only")
    d = pd.concat([native, hko], ignore_index=True)
    d["lon"] = d.lon % 360
    d = d[(d.lon >= 100) & (d.lon <= 180) & (d.lat >= 0) & (d.lat <= 40)]
    d = d.sort_values(["agency", "sid", "iso_time"])
    rows = []
    for (ag, sid), g in d.groupby(["agency", "sid"], sort=False):
        w = g.wind.to_numpy(float)
        if np.all(np.isnan(w)):
            continue
        idx = np.flatnonzero(w == np.nanmax(w))
        lat = g.lat.to_numpy(float)
        rows.append(dict(agency=ag, sid=sid, season=int(g.season.iloc[0]), genesis_lat=lat[0],
                         lmi_first=lat[idx[0]], lmi_last=lat[idx[-1]], lmi_mean=lat[idx].mean(),
                         n_tied=len(idx)))
    return pd.DataFrame(rows)

def annual(s, col):
    a = s.groupby("season")[["genesis_lat", col]].mean()
    a["disp"] = a[col] - a.genesis_lat
    return a

def run(src: Path, out: Path):
    s = storms(src)
    res = []
    for ag, g in s.groupby("agency"):
        for rule in ["lmi_first", "lmi_last", "lmi_mean"]:
            a = annual(g, rule)
            early, late = a.loc[1966:1995], a.loc[1996:2025]
            d = late.mean() - early.mean()
            res.append(dict(record=ag, tie_rule=rule, measure="two_period_difference",
                            n_years_early=len(early), n_years_late=len(late),
                            d_lmi_lat=d[rule], d_genesis_lat=d.genesis_lat, d_displacement=d.disp,
                            displacement_share=d.disp / d[rule]))
            for y0, y1 in [(1966, 2016), (1966, 2025)]:
                b = a.loc[y0:y1]
                t = b.index.to_numpy(float)
                sl = {c: np.polyfit(t, b[c].to_numpy(float), 1)[0] * 10 for c in [rule, "genesis_lat", "disp"]}
                res.append(dict(record=ag, tie_rule=rule, measure=f"trend_{y0}_{y1}_per_decade",
                                n_years_early=len(b), n_years_late=np.nan,
                                d_lmi_lat=sl[rule], d_genesis_lat=sl["genesis_lat"], d_displacement=sl["disp"],
                                displacement_share=sl["disp"] / sl[rule]))
    out_df = pd.DataFrame(res)
    cov = s.groupby("agency").agg(storms=("sid", "size"), first_season=("season", "min"), mean_tied=("n_tied", "mean"))
    out_df = out_df.merge(cov, left_on="record", right_index=True)
    out_df.to_csv(out, index=False, float_format="%.6f", lineterminator="\n", mode="x")
    print(out_df.round(3).to_string())

def main():
    parser = argparse.ArgumentParser(description="Split LMI latitude into genesis latitude and displacement.")
    parser.add_argument("tracks", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.tracks, args.output)


if __name__ == "__main__":
    main()
