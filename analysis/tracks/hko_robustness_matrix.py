"""HKO annual track-composition robustness configurations on the 2.5-degree grid.

Uses annual normalization, track-point and storm-equal weighting, endpoint and
leave-one-decade-out choices, signed projection, and 1,999 three-year block
permutations with seed 202406. Existing non-HKO CSV records are preserved.
"""
from pathlib import Path
import argparse,csv,io
import numpy as np
import pandas as pd
from core.paths import PROJECT
from analysis.shared import agency_catalogues
START,END,REGION,SEED,NPERM=1966,2025,(100.,180.,0.,40.),202406,1999
ENDPOINTS=[('end_2024_drop_1995',np.arange(1966,2025),(1966,1994),(1996,2024)),('start_1967_drop_1996',np.arange(1967,2026),(1967,1995),(1997,2025)),('exclude_2020_2025',np.arange(1966,2020),(1966,1992),(1993,2019))]
REMOVALS={'drop_1970s':(1970,1979),'drop_1980s':(1980,1989),'drop_1990s':(1990,1999),'drop_2000s':(2000,2009),'drop_2010s':(2010,2019),'drop_2020_2025':(2020,2025)}

def grid_edges(width):
    lon0, lon1, lat0, lat1 = REGION
    return np.arange(lon0, lon1 + width, width), np.arange(lat0, lat1 + width, width)

def annual_fields(tracks, years, width=2.5):
    lon_edges, lat_edges = grid_edges(width)
    shape = (len(lat_edges) - 1, len(lon_edges) - 1)
    point_fields, storm_fields = [], []
    for year in years:
        annual = tracks.loc[tracks["season"] == year]
        if annual.empty:
            raise ValueError(f"no eligible tracks in {year}")
        point = np.histogram2d(annual["lon"], annual["lat"], bins=[lon_edges, lat_edges])[0].T
        if point.sum() == 0:
            raise ValueError(f"no in-domain points in {year}")
        point_fields.append((point / point.sum()).reshape(-1))
        storm = np.zeros(shape, dtype=float)
        for _, group in annual.groupby("sid"):
            one = np.histogram2d(group["lon"], group["lat"], bins=[lon_edges, lat_edges])[0].T
            if one.sum() > 0:
                storm += one / one.sum()
        if storm.sum() == 0:
            raise ValueError(f"no valid storms in {year}")
        storm_fields.append((storm / storm.sum()).reshape(-1))
    return np.asarray(point_fields), np.asarray(storm_fields)

def temporal_blocks(length, block, years=None):
    positions = np.arange(length)
    if years is None:
        runs = [positions]
    else:
        years = np.asarray(years, dtype=int)
        runs = np.split(positions, np.flatnonzero(np.diff(years) != 1) + 1)
    return [run[i: i + block] for run in runs for i in range(0, len(run), block)]

def global_tv_permutation(fields, early_idx, late_idx, block, nperm, seed=SEED, years=None):
    change = fields[late_idx].mean(axis=0) - fields[early_idx].mean(axis=0)
    observed = float(0.5 * np.abs(change).sum())
    blocks = temporal_blocks(len(fields), block, years)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(nperm):
        order = np.concatenate([blocks[i] for i in rng.permutation(len(blocks))])
        diff = fields[order[late_idx]].mean(axis=0) - fields[order[early_idx]].mean(axis=0)
        exceed += 0.5 * np.abs(diff).sum() >= observed
    return observed, (exceed + 1) / (nperm + 1), change

def idx(years, early, late):
    return (np.flatnonzero((years >= early[0]) & (years <= early[1])),
            np.flatnonzero((years >= late[0]) & (years <= late[1])))

def rows_for(catalog, tracks, d_primary):
    base = np.arange(START, END + 1)
    point, storm = annual_fields(tracks, base)
    full = {"track_point": point, "storm_equal": storm}
    out = []
    def test(analysis, weighting, years, fields, early, late):
        e, l = idx(years, early, late)
        tv, p, change = global_tv_permutation(fields, e, l, block=3, nperm=NPERM, years=years)
        out.append(dict(analysis=analysis, catalog=catalog, weighting=weighting, grid_deg=2.5, block_years=3,
                        early_start=early[0], early_end=early[1], late_start=late[0], late_end=late[1],
                        n_early=len(e), n_late=len(l), tv=tv, block_permutation_p=p, n_permutations=NPERM,
                        signed_projection=float(change @ d_primary) / float(d_primary @ d_primary),
                        spatial_r_vs_primary=float(np.corrcoef(change, d_primary)[0, 1])))
    for w, f in full.items():
        test("cross_agency_primary", w, base, f, (1966, 1995), (1996, 2025))
    for label, years, early, late in ENDPOINTS:
        keep = np.isin(base, years)
        for w, f in full.items():
            test(label, w, years, f[keep], early, late)
    for label, (a, b) in REMOVALS.items():
        keep = ~((base >= a) & (base <= b))
        for w, f in full.items():
            test(label, w, base[keep], f[keep], (1966, 1995), (1996, 2025))
    return full, out


def combine_rows(matrix_path,hko,output):
    raw=matrix_path.read_bytes();lines=raw.splitlines(keepends=True)
    records=list(csv.DictReader(io.StringIO(raw.decode('utf-8'))))
    keep=[line for line,row in zip(lines[1:],records) if row['catalog']!='HKO']
    newline='\r\n' if lines[0].endswith(b'\r\n') else '\n'
    head=b''.join([lines[0],*keep])
    if not head.endswith(b'\n'):head+=newline.encode()
    payload=head+hko.to_csv(index=False,header=False,lineterminator=newline).encode('utf-8')
    frame=pd.read_csv(io.BytesIO(payload))
    identity=['analysis','catalog','weighting','grid_deg','block_years','early_start','early_end','late_start','late_end']
    if frame.duplicated(identity).any():raise ValueError('Duplicate robustness configuration')
    output.write_bytes(payload)
    return frame


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('data_root',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();root=a.data_root.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    matrix=root/'results/tracks/robustness_matrix_signed.csv';stored=pd.read_csv(matrix)
    with np.load(root/'data/reference/environment_validation/supplemental_validation/results/wnp_tc_redistribution_pattern.npz',allow_pickle=True) as z:
        years=z['years'];fields=z['PRIMARY_track_point_annual_fields']
        reference=fields[years>=1996].mean(0)-fields[years<=1995].mean(0)
    agency_catalogues.PROJECT=root
    agency_catalogues.output_path=lambda relative:root/relative
    tracks=agency_catalogues.all_catalogues()['HKO']['ts_points']
    tracks=tracks.loc[tracks.season.between(START,END)&tracks.lon.between(100,180)&tracks.lat.between(0,40)].sort_values(['season','sid','iso_time'])
    _,rows=rows_for('HKO',tracks,reference);hko=pd.DataFrame(rows)
    hko['spatial_corr_with_primary_same_grid']=np.nan;hko['projection_status']='computed'
    for col in ['signed_projection','spatial_r_vs_primary']:hko[col]=hko[col].round(4)
    hko=hko[stored.columns];hko.to_csv(a.output/'hko_robustness_rows.csv',index=False)
    combine_rows(matrix,hko,a.output/'robustness_matrix_signed.csv')

if __name__=='__main__':main()
