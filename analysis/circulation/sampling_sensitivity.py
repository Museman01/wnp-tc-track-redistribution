"""块置换零假设检验与分箱／起始年敏感性（配合 storm_time_sampling.py）。

自助给的是估计量的抽样区间，不是零假设检验。这里补一个与正文同构的
**3 年块置换**：把 60 个年份的时段标签按 3 年块打乱，重算各分量，得到双侧 p。
另外重跑 5.0° 分箱与 1982 起始年两组敏感性，检查符号是否稳定。
"""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import importlib.util
import json
import argparse
import csv
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np

HERE = _LAYOUT_ROOT / 'analysis/circulation'
TASK = _LAYOUT_ROOT

spec = importlib.util.spec_from_file_location("sts", _LAYOUT_ROOT / 'analysis/circulation/storm_time_sampling.py')
sts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sts)

NPERM = 9999
KEYS = ("total_change", "field", "sampling", "cosampling")
FIXED_KEYS = ("field", "seasonal", "cosampling")


def run_case(cfg, bin_deg, start_year, months):
    cfg = json.loads(json.dumps(cfg))
    cfg["grids"]["density_bin"] = bin_deg
    years = list(range(start_year, 2026))
    split = start_year + (2026 - start_year) // 2 - 1
    yp1 = [y for y in years if y <= split]
    yp2 = [y for y in years if y > split]

    grid = sts.coarsen_steering(cfg, months, years)
    U = sts.field_by_year_month(grid, months, years)
    valid = np.isfinite(U["u"]).all(axis=0) & np.isfinite(U["v"]).all(axis=0)
    Uf = {c: np.nan_to_num(U[c]) for c in ("u", "v")}

    frame = sts.read_ibtracs_agencies(
        sts.PROJECT / "data" / "raw" / "ibtracs_wp.csv", start=start_year, end=2025
    )
    case = {"bin_deg": bin_deg, "period_split": [start_year, split, split + 1, 2025],
            "agencies": {}}
    for agency in sts.AGENCIES:
        cat = sts.build_agency_catalog(frame, agency)
        dists, kept, _ = sts.annual_point_distributions(
            cat["tracks"], cfg, months, years, grid
        )
        a1 = [y for y in yp1 if y in dists]
        a2 = [y for y in yp2 if y in dists]
        rec = {}
        for comp in ("u", "v"):
            d = sts.decompose(Uf[comp], dists, years, a1, a2, valid)
            rec[comp] = {
                "total_change": d["total_change"],
                "field": d["field"],
                "sampling": d["sampling"],
                "cosampling": d["cosampling"],
                "fixed_space": d["fixed_space"],
            }
        case["agencies"][sts.AGENCY_LABEL[agency]] = rec
    return case


def permutation_test(cfg, months):
    years = list(range(1966, 2026))
    grid = sts.coarsen_steering_cached(cfg, months, years)
    U = sts.field_by_year_month(grid, months, years)
    valid = np.isfinite(U["u"]).all(axis=0) & np.isfinite(U["v"]).all(axis=0)
    Uf = {c: np.nan_to_num(U[c]) for c in ("u", "v")}
    frame = sts.read_ibtracs_agencies(
        sts.PROJECT / "data" / "raw" / "ibtracs_wp.csv", start=1966, end=2025
    )
    block = cfg["statistics"]["bootstrap_block"]
    rng = np.random.default_rng(cfg["statistics"]["random_seed"])

    blocks = [years[i:i + block] for i in range(0, len(years), block)]
    half = len(blocks) // 2

    out = {"n_permutations": NPERM, "block_years": block, "agencies": {}}
    for agency in sts.AGENCIES:
        cat = sts.build_agency_catalog(frame, agency)
        dists, kept, _ = sts.annual_point_distributions(
            cat["tracks"], cfg, months, years, grid
        )
        yp1 = [y for y in kept if y <= 1995]
        yp2 = [y for y in kept if y >= 1996]
        rec = {}
        for comp in ("u", "v"):
            obs = sts.decompose(Uf[comp], dists, years, yp1, yp2, valid)
            stats = {k: [] for k in KEYS}
            stats.update({f"fixed_{k}": [] for k in FIXED_KEYS})
            for _ in range(NPERM):
                order = rng.permutation(len(blocks))
                g1 = [y for b in order[:half] for y in blocks[b] if y in dists]
                g2 = [y for b in order[half:] for y in blocks[b] if y in dists]
                d = sts.decompose(Uf[comp], dists, years, g1, g2, valid)
                for k in KEYS:
                    stats[k].append(d[k])
                for k in FIXED_KEYS:
                    stats[f"fixed_{k}"].append(d["fixed_space"][k])
            res = {}
            for k in KEYS:
                null = np.asarray(stats[k])
                exceed = int(np.count_nonzero(np.abs(null) >= abs(obs[k])))
                res[k] = {
                    "observed": obs[k],
                    "p_two_sided": float((exceed + 1) / (NPERM + 1)),
                }
            for k in FIXED_KEYS:
                null = np.asarray(stats[f"fixed_{k}"])
                exceed = int(
                    np.count_nonzero(
                        np.abs(null) >= abs(obs["fixed_space"][k])
                    )
                )
                res[f"fixed_{k}"] = {
                    "observed": obs["fixed_space"][k],
                    "p_two_sided": float((exceed + 1) / (NPERM + 1)),
                }
            rec[comp] = res
        out["agencies"][sts.AGENCY_LABEL[agency]] = rec
    return out


def format_s7_pvalue(value):
    """Table S7: threshold at 0.001, otherwise decimal half-up to three places."""
    p = Decimal(str(value))
    return '<0.001' if p < Decimal('.001') else str(p.quantize(Decimal('.001'), rounding=ROUND_HALF_UP))


def table_s7_reporting(source_root):
    """Recover retained exceedance counts; return reporting rows and alternate-count differences."""
    root = Path(source_root)
    original = json.loads((root/'data/reference/storm_sampling/permutation_test.json').read_text())
    current = json.loads((root/'results/circulation/storm_sampling/permutation_test.json').read_text())
    if original['n_permutations'] != NPERM or current['n_permutations'] != NPERM:
        raise ValueError('Table S7 requires 9999 retained permutations')
    if original['block_years'] != 3 or current['block_years'] != 3:
        raise ValueError('Table S7 requires three-year blocks')
    with (root/'results/agency_comparison/table_s7_permutation_four_agency.csv').open(newline='') as f:
        vector = list(csv.DictReader(f))
    expected = {(a,c,k) for a in ['USA','JMA','CMA','HKO'] for c in ['u','v']
                for k in list(KEYS)+['fixed_field','fixed_seasonal','fixed_cosampling']}
    rows, differences, seen = [], [], set()
    for item in vector:
        agency = 'JMA' if item['agency'] == 'TOKYO' else item['agency']
        comp, term = item['wind_component'], item['term']
        key = (agency,comp,term)
        if key not in expected or key in seen:
            raise ValueError(f'Unexpected or repeated Table S7 term: {key}')
        seen.add(key)
        pv = float(item['p'])
        vector_count = round(pv*(NPERM+1))-1
        if abs(pv*(NPERM+1)-(vector_count+1)) > 1e-8:
            raise ValueError(f'Nonintegral HKO/vectorized count: {key}')
        if agency == 'HKO':
            count, source = vector_count, 'existing_HKO_count'
        else:
            ref = original['agencies'][agency][comp][term]
            raw_count = ref['p_two_sided']*NPERM
            count, source = round(raw_count), 'original_three_agency_count'
            if abs(raw_count-count) > 1e-8 or abs(float(item['change_ms'])-ref['observed']) > 1e-9:
                raise ValueError(f'Invalid retained count or effect: {key}')
            if abs(count-vector_count) > 1:
                raise ValueError(f'Exceedance counts differ by more than one: {key}')
            if (count+1)/(NPERM+1) != current['agencies'][agency][comp][term]['p_two_sided']:
                raise ValueError(f'Current add-one result disagrees with retained count: {key}')
        if not 0 <= count <= NPERM:
            raise ValueError(f'Count out of range: {key}')
        p = (count+1)/(NPERM+1)
        row = dict(agency=agency,wind_component=comp,term=term,change_ms=item['change_ms'],
                   n_permutations=NPERM,n_extreme=count,p_add_one=f'{p:.4f}',
                   count_source=source,p_display=format_s7_pvalue(f'{p:.4f}'))
        rows.append(row)
        if count != vector_count:
            differences.append(dict(agency=agency,component=comp,term=term,retained_count=count,
                                    vector_count=vector_count,retained_p=f'{p:.4f}',vector_p=item['p']))
    if seen != expected:
        raise ValueError('Missing Table S7 terms')
    return rows, differences


def write_s7_reporting(source_root, output_csv, overwrite=False):
    """Write Table S7 from exceedance counts using decimal half-up rounding."""
    destination = Path(output_csv)
    rows, differences = table_s7_reporting(source_root)
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('w' if overwrite else 'x',encoding='utf-8',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows, differences


def run_analysis(output_dir):
    dest = (_LAYOUT_ROOT / output_dir).resolve()
    if not dest.is_relative_to(_LAYOUT_ROOT.resolve()): raise ValueError('Output must remain within the project')
    if dest.exists(): raise FileExistsError('Provide a new output directory')
    dest.mkdir(parents=True)
    cfg = sts.load_config()
    months = cfg["typhoon_season"]
    perm = permutation_test(cfg, months)
    with open(dest / "permutation_test.json", "w", encoding="utf-8") as fh:
        json.dump(perm, fh, indent=2, ensure_ascii=False)

    sens = {
        "bin_5deg_1966": run_case(cfg, 5.0, 1966, months),
        "bin_2p5deg_1982": run_case(cfg, 2.5, 1982, months),
    }
    with open(dest / "sensitivity.json", "w", encoding="utf-8") as fh:
        json.dump(sens, fh, indent=2, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-s7',action='store_true',help='Write Table S7 from retained exceedance counts')
    parser.add_argument('--source-root',type=Path,default=TASK)
    parser.add_argument('--output-csv',type=Path)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--output-dir',type=Path)
    args = parser.parse_args()
    if args.report_s7:
        if not args.execute or args.output_csv is None:
            parser.error('--report-s7 requires --execute and an explicit --output-csv')
        if args.output_dir is not None: parser.error('--output-dir applies to the sampling analysis')
        output=(_LAYOUT_ROOT / args.output_csv).resolve()
        if not output.is_relative_to(_LAYOUT_ROOT.resolve()): parser.error('output must remain within the project')
        rows, differences = write_s7_reporting(args.source_root,output)
        print(str(args.output_csv))
    else:
        if args.output_csv is not None or not args.execute or args.output_dir is None or args.source_root != TASK:
            parser.error('sampling analysis requires --execute and --output-dir')
        run_analysis(args.output_dir)


if __name__ == '__main__':
    main()
