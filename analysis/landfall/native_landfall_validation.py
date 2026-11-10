from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT



import argparse


from pathlib import Path



import geopandas as gpd


import numpy as np

import pandas as pd

import shapely

from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPoint, Point




from analysis.shared.agency_paths import PROJECT, load_config




from core.landfall_gshhg import attribute_coast

from analysis.shared.agency_data import AGENCIES, agency_flag, build_agency_catalog, native_ts_mask, read_ibtracs_agencies



from analysis.shared.block_resampling import bh_fdr, block_bootstrap_many, block_permutation_many, trend_summary





ANALYSIS = _LAYOUT_ROOT / 'results/landfall/native_validation'

AUTHORITY_PATH = (

    _LAYOUT_ROOT / 'data/processed/exclusive_coast/classified_exact_vector_events_admin0_corrected.csv'

)







class LocalPolygonSet:

    """Small vectorized land predicate without constructing polygon unions."""



    def __init__(self, geometries):

        self.geometries = geometries



    def covers(self, point):

        return bool(np.any(shapely.covers(self.geometries, point)))





def geometry_points(geometry):

    """Return point-like boundary intersections, including overlap endpoints."""

    if geometry is None or geometry.is_empty:

        return []

    if isinstance(geometry, Point):

        return [geometry]

    if isinstance(geometry, (MultiPoint, GeometryCollection)):

        out = []

        for part in geometry.geoms:

            out.extend(geometry_points(part))

        return out

    if isinstance(geometry, LineString):

        coords = list(geometry.coords)

        return [Point(coords[0]), Point(coords[-1])] if coords else []

    if isinstance(geometry, MultiLineString):

        out = []

        for part in geometry.geoms:

            out.extend(geometry_points(part))

        return out

    return []





def sea_to_land_fractions(line: LineString, intersection, prepared_land) -> list[float]:

    fractions = sorted(line.project(p, normalized=True) for p in geometry_points(intersection))

    unique = []

    for value in fractions:

        value = float(np.clip(value, 0.0, 1.0))

        if not unique or abs(value - unique[-1]) > 1e-8:

            unique.append(value)

    if not unique:

        return []

    bounds = [0.0, *unique, 1.0]

    states = []

    for left, right in zip(bounds[:-1], bounds[1:]):

        fraction = (left + right) / 2 if right - left > 1e-10 else min(1.0, left + 1e-8)

        states.append(bool(prepared_land.covers(line.interpolate(fraction, normalized=True))))

    return [fraction for i, fraction in enumerate(unique) if not states[i] and states[i + 1]]





def build_segments(tracks: pd.DataFrame, max_gap_hours: float) -> pd.DataFrame:

    frame = tracks.sort_values(["sid", "iso_time"]).copy()

    for column in ["sid", "iso_time", "lat", "lon", "wind"]:

        frame[f"next_{column}"] = frame.groupby("sid")[column].shift(-1)

    frame["gap_hours"] = (frame.next_iso_time - frame.iso_time).dt.total_seconds() / 3600

    return frame.loc[

        frame.sid.eq(frame.next_sid)

        & frame.gap_hours.gt(0)

        & frame.gap_hours.le(max_gap_hours)

        & frame[["lat", "lon", "next_lat", "next_lon"]].notna().all(axis=1)

    ].copy()





def refine_candidate_crossings(tracks, agency, candidates, coast_geometries, coast_tree, cfg):

    """Refine conservative raster candidates to exact vector-coast crossings.



    The 0.02-degree mask is used only to identify a small set of candidate

    segments. Final position, time and wind are recomputed from the exact

    intersection with the GSHHG high-resolution polygon boundary.

    """

    segments = build_segments(tracks, cfg["landfall"]["maximum_track_gap_hours"])

    by_sid = {sid: group for sid, group in segments.groupby("sid", sort=False)}

    events, unresolved, cache = [], [], {}

    for candidate in candidates.sort_values("time").itertuples():

        group = by_sid.get(candidate.sid)

        if group is None:

            unresolved.append({"agency": agency, "sid": candidate.sid, "time": candidate.time, "reason": "sid_not_in_track_segments"})

            continue

        matched = group.loc[group.iso_time.le(candidate.time) & group.next_iso_time.ge(candidate.time)]

        if matched.empty:

            unresolved.append({"agency": agency, "sid": candidate.sid, "time": candidate.time, "reason": "candidate_not_bracketed"})

            continue

        row = matched.iloc[0]

        cache_key = (candidate.sid, row.iso_time, row.next_iso_time)

        if cache_key not in cache:

            line = LineString([(row.lon, row.lat), (row.next_lon, row.next_lat)])

            polygon_indices = coast_tree.query(line)

            if len(polygon_indices):

                local_geometries = coast_geometries[polygon_indices]

                intersections = shapely.intersection(shapely.boundary(local_geometries), line)

                intersection = GeometryCollection([item for item in intersections if not item.is_empty])

                fractions = sea_to_land_fractions(line, intersection, LocalPolygonSet(local_geometries))

            else:

                fractions = []

            cache[cache_key] = (line, fractions)

        line, fractions = cache[cache_key]

        if not fractions:

            unresolved.append({"agency": agency, "sid": candidate.sid, "time": candidate.time, "reason": "no_exact_sea_to_land_intersection"})

            continue

        approximate_fraction = (candidate.time - row.iso_time) / (row.next_iso_time - row.iso_time)

        fraction = min(fractions, key=lambda value: abs(value - approximate_fraction))

        event_time = row.iso_time + (row.next_iso_time - row.iso_time) * fraction

        lat = row.lat + fraction * (row.next_lat - row.lat)

        lon = row.lon + fraction * (row.next_lon - row.lon)

        wind = np.nan

        if pd.notna(row.wind) and pd.notna(row.next_wind):

            wind = row.wind + fraction * (row.next_wind - row.wind)

        events.append({

            "agency": agency, "sid": row.sid, "season": int(row.season),

            "time": event_time, "lat": float(lat), "lon": float(lon), "wind": float(wind),

            "coast": attribute_coast(float(lat), float(lon)), "segment_start": row.iso_time,

            "segment_end": row.next_iso_time, "start_lat": row.lat, "start_lon": row.lon,

            "end_lat": row.next_lat, "end_lon": row.next_lon, "segment_fraction": fraction,

            "candidate_source": "0.02-degree raster sea-to-land screen",

        })

    events = pd.DataFrame(events)

    if events.empty:

        return events, pd.DataFrame(unresolved)

    events = events.sort_values("time").drop_duplicates(["sid", "time"])

    kept = []

    minimum = cfg["landfall"]["minimum_event_separation_hours"]

    for sid, group in events.sort_values("time").groupby("sid", sort=False):

        last = None

        for index, row in group.iterrows():

            if last is None or (row.time - last).total_seconds() / 3600 >= minimum:

                kept.append(index)

                last = row.time

    events = events.loc[kept].sort_values(["season", "sid", "time"]).reset_index(drop=True)

    return events, pd.DataFrame(unresolved)





def definitions(events: pd.DataFrame) -> dict[str, pd.DataFrame]:

    ordered = events.sort_values(["sid", "time"])

    return {

        "first_landfall": ordered.drop_duplicates("sid", keep="first"),

        "all_events": ordered,

        "unique_storm_coast": ordered.drop_duplicates(["sid", "coast"], keep="first"),

    }





def attach_native_stage_and_exclusive_coast(

    events: pd.DataFrame,

    source: pd.DataFrame,

) -> pd.DataFrame:

    """Attach agency-native TS state and the current exclusive coast authority."""

    authority_path = AUTHORITY_PATH

    authority = pd.read_csv(

        authority_path,

        usecols=[

            "agency",

            "sid",

            "segment_start",

            "segment_end",

            "coast_exclusive",

            "nearest_segment",

            "nearest_distance_km",

        ],

        parse_dates=["segment_start", "segment_end"],

    )

    keys = ["agency", "sid", "segment_start", "segment_end"]

    if authority.duplicated(keys).any():

        raise RuntimeError("Exclusive-coast authority contains duplicate event keys")



    out = events.copy()

    out["segment_start"] = pd.to_datetime(out["segment_start"])

    out["segment_end"] = pd.to_datetime(out["segment_end"])

    out = out.merge(authority, on=keys, how="left", validate="one_to_one")

    if out["coast_exclusive"].isna().any():

        raise RuntimeError("Exact events did not fully match exclusive-coast authority")



    chunks = []

    for agency, info in AGENCIES.items():

        positioned = (

            agency_flag(source, agency, intensity=False)

            & source[info["lat"]].notna()

            & source[info["lon"]].notna()

        )

        state = source.loc[positioned, ["SID", "ISO_TIME"]].copy()

        state.columns = ["sid", "state_time"]

        state["native_ts"] = native_ts_mask(source, agency).loc[positioned].to_numpy(bool)

        state["intensity_original"] = (

            agency_flag(source, agency, intensity=True)

            .loc[positioned]

            .to_numpy(bool)

        )

        state = state.sort_values(["sid", "state_time"]).drop_duplicates(

            ["sid", "state_time"], keep="last"

        )

        group = out.loc[out["agency"].eq(agency)].copy()

        for side, time_column in (("start", "segment_start"), ("end", "segment_end")):

            renamed = state.rename(

                columns={

                    "state_time": time_column,

                    "native_ts": f"{side}_native_ts",

                    "intensity_original": f"{side}_intensity_original",

                }

            )

            group = group.merge(

                renamed,

                on=["sid", time_column],

                how="left",

                validate="many_to_one",

            )

        chunks.append(group)

    out = pd.concat(chunks, ignore_index=True)

    endpoint_columns = [

        "start_native_ts",

        "end_native_ts",

        "start_intensity_original",

        "end_intensity_original",

    ]

    if out[endpoint_columns].isna().any().any():

        raise RuntimeError("At least one exact event lacks an agency-native endpoint state")

    if not (

        out["start_intensity_original"].astype(bool)

        & out["end_intensity_original"].astype(bool)

    ).all():

        raise RuntimeError("At least one endpoint is not an original/verified intensity report")

    out["pre_crossing_native_ts"] = out["start_native_ts"].astype(bool)

    out["either_endpoint_native_ts"] = (

        out["start_native_ts"].astype(bool) | out["end_native_ts"].astype(bool)

    )

    out["both_endpoints_native_ts"] = (

        out["start_native_ts"].astype(bool) & out["end_native_ts"].astype(bool)

    )

    out["coast_legacy"] = out["coast"]

    out["coast"] = out["coast_exclusive"]

    return out.sort_values(["agency", "season", "sid", "time"]).reset_index(drop=True)





def annual_metrics(events: pd.DataFrame, years: np.ndarray) -> dict[str, float]:

    rows = []

    for year in years:

        values = events.loc[events.season.eq(year), "lat"].to_numpy(float)

        rows.append({

            "year": year, "n_events": len(values), "n_storms": events.loc[events.season.eq(year), "sid"].nunique(),

            "mean_lat": np.mean(values) if len(values) else np.nan,

            "median_lat": np.median(values) if len(values) else np.nan,

            "q25_lat": np.quantile(values, 0.25) if len(values) else np.nan,

            "q75_lat": np.quantile(values, 0.75) if len(values) else np.nan,

        })

    return rows





def summarize(annual: pd.DataFrame, cfg: dict) -> pd.DataFrame:

    rows = []

    metrics = ["mean_lat", "median_lat", "q25_lat", "q75_lat"]

    for (agency, definition), group in annual.groupby(["agency", "definition"]):

        group = group.sort_values("year")

        for start, end in [(1966, 2025), (1966, 2024), (1982, 2025), (1982, 2024)]:

            period = group.loc[group.year.between(start, end)]

            for metric in metrics:

                trend = trend_summary(period.year, period[metric])

                row = {"agency": agency, "definition": definition, "metric": metric, "start": start, "end": end, **trend,

                       "period_difference": np.nan, "period_ci_low": np.nan, "period_ci_high": np.nan, "period_block_p": np.nan}

                rows.append(row)

    out = pd.DataFrame(rows)

    for end in [2024, 2025]:

        target = out.loc[(out.start == 1966) & (out.end == end)]

        series = []

        for r in target.itertuples():

            d = annual.loc[annual.agency.eq(r.agency) & annual.definition.eq(r.definition) & annual.year.between(1966, end)].sort_values("year")

            series.append(d[r.metric].to_numpy(float))

        matrix = np.column_stack(series)

        years = np.arange(1966, end + 1)

        early_idx = np.flatnonzero(years <= 1995)

        late_idx = np.flatnonzero(years >= 1996)

        difference, p = block_permutation_many(matrix, early_idx, late_idx, block=3, nperm=cfg["n_permutations"], seed=cfg["random_seed"])

        lo, hi = block_bootstrap_many(matrix[early_idx], matrix[late_idx], block=3, nboot=cfg["n_bootstrap"], seed=cfg["random_seed"])

        out.loc[target.index, "period_difference"] = difference

        out.loc[target.index, "period_ci_low"] = lo

        out.loc[target.index, "period_ci_high"] = hi

        out.loc[target.index, "period_block_p"] = p

    out["mk_fdr_family"] = (

        "landfall_latitude_trend_"

        + out.definition

        + "_"

        + out.metric

        + "_"

        + out.start.astype(str)

        + "_"

        + out.end.astype(str)

    )

    out["mk_q_bh"] = np.nan

    for _, idx in out.groupby("mk_fdr_family").groups.items():

        out.loc[idx, "mk_q_bh"] = bh_fdr(out.loc[idx, "mk_p"])

    out["period_fdr_family"] = np.where(

        out.period_block_p.notna(),

        "landfall_latitude_period_"

        + out.definition

        + "_"

        + out.metric

        + "_1966_"

        + out.end.astype(str),

        "",

    )

    out["period_q_bh"] = np.nan

    for family, idx in out.loc[out.period_block_p.notna()].groupby("period_fdr_family").groups.items():

        out.loc[idx, "period_q_bh"] = bh_fdr(out.loc[idx, "period_block_p"])

    return out





def leave_one_coast_out(events: pd.DataFrame, cfg: dict) -> pd.DataFrame:

    years = np.arange(1966, 2026)

    rows, series = [], []

    for agency, agency_events in events.groupby("agency"):

        for definition, defined in definitions(agency_events).items():

            if definition not in {"first_landfall", "all_events"}:

                continue

            for excluded in ["none", *sorted(defined.coast.unique())]:

                sample = defined if excluded == "none" else defined.loc[defined.coast.ne(excluded)]

                annual = pd.DataFrame(annual_metrics(sample, years))

                early = np.flatnonzero(annual.year.between(1966, 1995))

                late = np.flatnonzero(annual.year.between(1996, 2025))

                rows.append({"agency": agency, "definition": definition, "excluded_coast": excluded, "n_events": len(sample)})

                series.append(annual.mean_lat.to_numpy(float))

    matrix = np.column_stack(series)

    difference, p = block_permutation_many(matrix, early, late, block=3, nperm=cfg["n_permutations"], seed=cfg["random_seed"])

    out = pd.DataFrame(rows)

    out["period_difference_mean_lat"] = difference

    out["block_permutation_p"] = p

    return out





def run(*, reuse_existing: bool = False, candidate_events: Path | None = None, output_dir: Path = ANALYSIS):
    if not reuse_existing and candidate_events is None:
        raise ValueError('Provide candidate events or use the existing exact crossings')
    out = (PROJECT / output_dir).resolve()
    if not out.is_relative_to(PROJECT.resolve()):
        raise ValueError('Output directory must remain within the project')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'Refusing to overwrite existing results: {out}')
    out.mkdir(parents=True, exist_ok=True)

    cfg = load_config()

    source = read_ibtracs_agencies(PROJECT / "data/raw/IBTrACS.WP.v04r01.csv", start=1966, end=2025)

    full_lifecycle_path = out / "landfall_events_exact_full_lifecycle.csv"

    if reuse_existing:
        source_path = PROJECT / 'results/landfall/landfall_events_exact_full_lifecycle.csv'
        full_lifecycle = pd.read_csv(source_path, parse_dates=["time", "segment_start", "segment_end"])

    else:

        catalogs = {agency: build_agency_catalog(source, agency) for agency in AGENCIES}

        candidate_events = pd.read_csv(candidate_events, parse_dates=["time"])

        candidate_events = candidate_events.loc[candidate_events.season.between(1966, 2025)]

        coast = gpd.read_file(cfg["landfall"]["coastline"], bbox=(80, 0, 180, 65))

        coast_geometries = coast.geometry.to_numpy()

        coast_tree = shapely.STRtree(coast_geometries)

        all_events, all_unresolved = [], []

        for agency, catalog in catalogs.items():

            events_agency, unresolved_agency = refine_candidate_crossings(

                catalog["tracks"], agency,

                candidate_events.loc[candidate_events.agency.eq(agency)],

                coast_geometries, coast_tree, cfg,

            )

            all_events.append(events_agency)

            all_unresolved.append(unresolved_agency)

        full_lifecycle = pd.concat(all_events, ignore_index=True)

        unresolved = pd.concat(all_unresolved, ignore_index=True)

        unresolved.to_csv(out / "unresolved_candidates.csv", index=False)



    if not reuse_existing:
        full_lifecycle.to_csv(full_lifecycle_path, index=False)

    stage_audit = attach_native_stage_and_exclusive_coast(full_lifecycle, source)

    stage_audit.to_csv(out / "landfall_events_exact_stage_audit.csv", index=False)

    events = stage_audit.loc[stage_audit["pre_crossing_native_ts"]].copy()

    events.to_csv(out / "landfall_events_exact.csv", index=False)

    years = np.arange(1966, 2026)

    annual_rows = []

    for agency, group in events.groupby("agency"):

        for definition, sample in definitions(group).items():

            for row in annual_metrics(sample, years):

                annual_rows.append({"agency": agency, "definition": definition, **row})

    annual = pd.DataFrame(annual_rows)

    annual.to_csv(out / "landfall_latitude_annual.csv", index=False)

    summary = summarize(annual, cfg)

    summary.to_csv(out / "landfall_latitude_summary.csv", index=False)

    leave = leave_one_coast_out(events, cfg)

    leave.to_csv(out / "leave_one_coast_out.csv", index=False)

if __name__ == "__main__":

    import sys

    if '--execute' not in sys.argv:

        raise SystemExit('Explicit execution required: add --execute. Output writing requires explicit execution.')

    sys.argv.remove('--execute')

    parser = argparse.ArgumentParser()

    parser.add_argument(

        "--reuse-existing",

        action="store_true",

        help="Reuse the existing exact crossing table and rerun stage/statistical layers.",

    )

    parser.add_argument("--candidate-events", type=Path)
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS)
    args = parser.parse_args()
    run(reuse_existing=args.reuse_existing, candidate_events=args.candidate_events, output_dir=args.output_dir)

