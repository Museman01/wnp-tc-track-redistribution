"""Annual agency track composition, map correlation and numerical reference comparison."""
from __future__ import annotations
import numpy as np
import pandas as pd


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


def numeric_comparison(new, old, keys, columns, source):

    merged = new.merge(old, on=keys, suffixes=("_new", "_reference"), how="outer", indicator=True)

    rows = []

    for _, row in merged.iterrows():

        base = {key: row[key] for key in keys}

        for column in columns:

            left = row.get(f"{column}_new", np.nan)

            right = row.get(f"{column}_reference", np.nan)

            match = bool(

                row["_merge"] == "both"

                and ((pd.isna(left) and pd.isna(right)) or np.isclose(left, right, rtol=1e-10, atol=1e-12, equal_nan=True))

            )

            rows.append({**base, "source": source, "field": column, "new": left, "reference": right, "match": match})

    return rows
