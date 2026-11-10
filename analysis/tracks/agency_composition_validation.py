from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT



from itertools import combinations




import numpy as np

import pandas as pd



from analysis.shared.agency_paths import PROJECT, load_config



from core.utils import load_config as load_project_config

from analysis.shared.agency_data import AGENCIES, build_agency_catalog, read_ibtracs_agencies

from analysis.shared.inference_stats import add_family_fdr, compositional_change_test, trend_row





def annual_path_composition(points, years, lon_edges, lat_edges):

    points = points.loc[

        points["lon"].between(lon_edges[0], lon_edges[-1])

        & points["lat"].between(lat_edges[0], lat_edges[-1])

    ]

    rows = []

    for year in years:

        group = points.loc[points["season"].eq(year)]

        field = np.histogram2d(group["lon"], group["lat"], bins=[lon_edges, lat_edges])[0].T

        if field.sum() == 0:

            raise ValueError(f"no agency path points in {year}")

        rows.append((field / field.sum()).ravel())

    return np.asarray(rows)





def safe_corr(a, b):

    a, b = np.asarray(a, float), np.asarray(b, float)

    valid = np.isfinite(a) & np.isfinite(b)

    return float(np.corrcoef(a[valid], b[valid])[0, 1]) if valid.sum() >= 3 else np.nan





def run(output_dir):

    spec = load_config()

    project_cfg = load_project_config()

    from core.paths import create_new_output_dir
    out = create_new_output_dir(output_dir)

    years = np.arange(1966, 2026)

    lmi_years = np.arange(1982, 2026)

    raw = PROJECT / "data" / "raw" / "IBTrACS.WP.v04r01.csv"

    source = read_ibtracs_agencies(raw, start=1945, end=2025)

    catalogs = {agency: build_agency_catalog(source, agency) for agency in AGENCIES}

    common_lmi = set.intersection(*(set(catalogs[a]["lmi"]["sid"]) for a in AGENCIES))



    annual_rows = []

    trend_rows = []

    for agency, catalog in catalogs.items():

        frequency = catalog["frequency"].set_index("season")["n_tc"].reindex(years, fill_value=0)

        lmi = catalog["lmi"].loc[catalog["lmi"]["season"].between(1982, 2025)]

        full = lmi.groupby("season")["lmi_lat"].mean().reindex(lmi_years)

        common = lmi.loc[lmi["sid"].isin(common_lmi)].groupby("season")["lmi_lat"].mean().reindex(lmi_years)

        for year in years:

            annual_rows.append({

                "agency": agency, "season": year, "n_tc": int(frequency.loc[year]),

                "mean_lmi_lat_full": full.loc[year] if year >= 1982 else np.nan,

                "mean_lmi_lat_common": common.loc[year] if year >= 1982 else np.nan,

            })

        for end in (2024, 2025):

            for variable, start, metric, values in [

                ("n_tc", 1966, "frequency", frequency.loc[1966:end]),

                ("mean_lmi_lat_full", 1982, "lmi_full_catalog", full.loc[1982:end]),

                ("mean_lmi_lat_common", 1982, "lmi_common_storms", common.loc[1982:end]),

            ]:

                trend_rows.append(trend_row(

                    values.index, values.values, label=agency,

                    family=f"multiagency_{metric}_{start}_{end}", cfg=project_cfg,

                    extra={"analysis": metric, "agency": agency, "variable": variable},

                ))

    annual = pd.DataFrame(annual_rows)

    trends = add_family_fdr(trend_rows)

    annual.to_csv(out / "core_crossagency_annual.csv", index=False)

    trends.to_csv(out / "core_crossagency_recheck_trends.csv", index=False)



    lon_edges = np.arange(100, 180 + 2.5, 2.5)

    lat_edges = np.arange(0, 40 + 2.5, 2.5)

    compositions = {a: annual_path_composition(c["ts_points"], years, lon_edges, lat_edges) for a, c in catalogs.items()}

    results = {

        a: compositional_change_test(x, np.arange(30), np.arange(30, 60),

                                     nperm=spec["n_permutations"], block=3, seed=spec["random_seed"])

        for a, x in compositions.items()

    }

    redistribution = pd.DataFrame([{

        "agency": a, "analysis": "path_density", "early_start": 1966, "early_end": 1995,

        "late_start": 1996, "late_end": 2025, "total_variation": r["tv"],

        "block_permutation_p": r["global_p"], "n_permutations": r["nperm"], "block_years": r["block"],

    } for a, r in results.items()])

    redistribution.to_csv(out / "core_crossagency_recheck_redistribution.csv", index=False)



    agreement = pd.DataFrame([{

        "agency_left": left, "agency_right": right, "metric": "path_change_map",

        "correlation": safe_corr(results[left]["change"], results[right]["change"]),

    } for left, right in combinations(AGENCIES, 2)])

    agreement.to_csv(out / "core_crossagency_recheck_agreement.csv", index=False)

    np.savez_compressed(out / "annual_path_composition_2p5deg.npz", years=years, lon_edges=lon_edges,

                        lat_edges=lat_edges, **compositions)



if __name__ == "__main__":

    from pathlib import Path
    import argparse

    parser = argparse.ArgumentParser(description="Agency composition and core annual validation.")
    parser.add_argument("--execute", action="store_true", help="confirm execution")
    parser.add_argument("--output-dir", required=True, help="new project-local output directory")
    args = parser.parse_args()
    if not args.execute:
        parser.error("add --execute to write outputs")
    run(args.output_dir)
