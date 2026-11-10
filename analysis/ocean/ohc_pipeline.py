"""Run the OHC exposure analysis using configured scientific inputs."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from analysis.ocean.ohc_common import (
    DEFAULT_CONFIG,
    load_config,
    project_path,
    require_execute,
    resolve_config_path,
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--phase", choices=["prepare", "analyze"], required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--use-prepared-ocean", action="store_true", help="use the included regional ocean files in the prepare phase")
    args = parser.parse_args()
    config_path = resolve_config_path(args.config)
    config = load_config(config_path)
    require_execute(args.execute)
    stages = PREPARE_STAGES if args.phase == "prepare" else ANALYZE_STAGES
    if args.phase == "prepare" and args.use_prepared_ocean:
        prepared = project_path(config["ohc"]["prepared_dir"])
        years = range(int(config["study"]["start_year"]), int(config["study"]["end_year"]) + 1)
        missing = [year for year in years if not (prepared / f"oras5_ohc300_{year}.nc").is_file()]
        if missing:
            raise FileNotFoundError(f"regional ORAS5 inputs missing for years: {missing}")
        for key in ["ibtracs_csv", "primary_csv"]:
            if not project_path(config["tracks"][key]).is_file():
                raise FileNotFoundError(config["tracks"][key])
        stages = stages[1:]
    for stage in stages:
        command = [sys.executable, "-m", stage, "--config", str(config_path), "--execute"]
        if args.overwrite:
            command.append("--overwrite")
        subprocess.run(command, check=True, cwd=project_path("."))


if __name__ == "__main__":
    main()
