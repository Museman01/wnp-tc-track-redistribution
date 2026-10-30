"""Run the OHC exposure analysis using configured scientific inputs."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from analysis.ocean.ohc_common import (
    DEFAULT_CONFIG,
    PACKAGE_ROOT,
    load_config,
    project_path,
    require_execute,
    resolve_config_path,
    sha256_file,
    validate_product_seam,
)


PREPARE_STAGES = [
    'analysis.ocean.ohc_preparation',
    'analysis.ocean.ohc_tracks',
    'analysis.ocean.ohc_matching',
    'analysis.ocean.ohc_product_boundary',
]

ANALYZE_STAGES = [
    'analysis.ocean.ohc_decomposition',
]


def write_manifest(config: dict) -> None:
    results = _LAYOUT_ROOT / 'results/ocean/ohc300'
    manifest_path = project_path(config["pipeline"]["manifest_csv"])
    files = sorted(path for path in results.rglob("*") if path.is_file() and path != manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"])
        writer.writeheader()
        for path in files:
            writer.writerow({
                "path": str(path.relative_to(project_path("."))),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })


def check_inputs(config: dict, phase: str, use_prepared_ocean: bool) -> dict:
    """Read input schemas without running or overwriting scientific analyses."""
    if phase == "analyze":
        validate_product_seam(config)
        source = project_path(config["matching"]["output_csv"])
        import gzip
        opener = gzip.open if source.suffix == ".gz" else open
        with opener(source, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"agency", "definition", "sid", "season", "ohc_year", "ohc_month", "match_status", "ohc300", "point_is_land"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"matched-point input lacks columns: {sorted(missing)}")
            if next(reader, None) is None:
                raise ValueError("matched-point input is empty")
        return {"status": "PASS", "phase": phase, "matched_points": source.relative_to(project_path(".")).as_posix(), "product_boundary": "PASS"}
    if use_prepared_ocean:
        prepared = project_path(config["ohc"]["prepared_dir"])
        years = range(int(config["study"]["start_year"]), int(config["study"]["end_year"]) + 1)
        missing = [year for year in years if not (prepared / f"oras5_ohc300_{year}.nc").is_file()]
        if missing:
            raise FileNotFoundError(f"regional ORAS5 inputs missing for years: {missing}")
    else:
        for key in ["consolidated_dir", "operational_dir"]:
            folder = project_path(config["ohc"][key])
            if not folder.is_dir() or not any(folder.iterdir()):
                raise FileNotFoundError(f"global ORAS5 input directory is absent: {folder}; use --use-prepared-ocean for the included regional files")
    for key in ["ibtracs_csv", "primary_csv"]:
        if not project_path(config["tracks"][key]).is_file():
            raise FileNotFoundError(config["tracks"][key])
    return {"status": "PASS", "phase": phase, "regional_ocean_inputs": use_prepared_ocean}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--phase", choices=["prepare", "analyze"], required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--check-inputs", action="store_true", help="validate input availability and schemas without analysis")
    parser.add_argument("--use-prepared-ocean", action="store_true", help="use the included regional ocean files in the prepare phase")
    args = parser.parse_args()
    config_path = resolve_config_path(args.config)
    config = load_config(config_path)
    if args.check_inputs:
        print(json.dumps(check_inputs(config, args.phase, args.use_prepared_ocean), indent=2))
        return
    require_execute(args.execute)
    if args.phase == "analyze":
        validate_product_seam(config)

    stages = PREPARE_STAGES if args.phase == "prepare" else ANALYZE_STAGES
    if args.phase == "prepare" and args.use_prepared_ocean:
        check_inputs(config, args.phase, True)
        stages = stages[1:]
    for stage in stages:
        command = [sys.executable, str(_LAYOUT_ROOT / 'analysis/ocean' / stage), "--config", str(config_path), "--execute"]
        if args.overwrite:
            command.append("--overwrite")
        subprocess.run(command, check=True, cwd=project_path("."))
    write_manifest(config)


if __name__ == "__main__":
    main()
