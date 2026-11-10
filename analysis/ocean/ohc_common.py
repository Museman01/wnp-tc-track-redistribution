"""Scientific input validation and path handling for OHC exposure."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


PACKAGE_ROOT = _LAYOUT_ROOT
PROJECT_ROOT = _LAYOUT_ROOT
DEFAULT_CONFIG = _LAYOUT_ROOT / 'config/ohc.json'


def coordinate_grids_equivalent(
    current: np.ndarray,
    reference: np.ndarray,
    atol_degrees: float,
) -> bool:
    """Treat sub-grid float32 coordinate quantization as equivalent geometry."""
    current = np.asarray(current, dtype="float64")
    reference = np.asarray(reference, dtype="float64")
    if current.shape != reference.shape:
        return False
    return bool(
        np.allclose(
            current,
            reference,
            rtol=0.0,
            atol=float(atol_degrees),
            equal_nan=True,
        )
    )


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()






def load_config(path: Path | None = None) -> dict[str, Any]:
    config_path = resolve_config_path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != 1:
        raise ValueError("unsupported config schema_version")
    return config


def resolve_config_path(path: Path | None = None) -> Path:
    return Path(path or DEFAULT_CONFIG).resolve()


def project_path(value: str | Path) -> Path:
    candidate = Path(value)
    resolved = candidate.resolve() if candidate.is_absolute() else (PROJECT_ROOT / candidate).resolve()
    if resolved != PROJECT_ROOT and PROJECT_ROOT not in resolved.parents:
        raise ValueError(f"path escapes project root: {resolved}")
    return resolved


def require_execute(enabled: bool) -> None:
    if not enabled:
        raise SystemExit(
            "Explicit execution required: add --execute"
        )






def validate_product_seam(config: dict) -> dict:
    """Check numerical product-boundary diagnostics."""
    settings = config["seam_audit"]
    report_path = project_path(settings["output_json"])
    if not report_path.is_file():
        raise FileNotFoundError(f"product-boundary diagnostics are absent: {report_path}")
    with report_path.open("r", encoding="utf-8") as handle:
        report = json.load(handle)
    if report.get("structural_errors") != []:
        raise ValueError("ORAS5 product boundary has missing or failed structural diagnostics")
    left, right = int(settings["left_year"]), int(settings["right_year"])
    if right != left + 1 or (report.get("left_year"), report.get("right_year")) != (left, right):
        raise ValueError("product-boundary years do not match the scientific configuration")
    tolerance = float(config["ohc"]["coordinate_equality_atol_degrees"])
    if float(report["coordinate_tolerance_degrees"]) != tolerance:
        raise ValueError("product-boundary coordinate tolerance differs from the configuration")
    for key in ["seam_max_abs_latitude_difference_degrees", "seam_max_abs_longitude_difference_degrees"]:
        value = float(report[key])
        if not np.isfinite(value) or value < 0.0 or value > tolerance:
            raise ValueError(f"product-boundary coordinate mismatch: {key}")
    for key in ["seam_changed_finite_any_cells", "seam_changed_finite_all_cells"]:
        if report.get(key) != 0:
            raise ValueError(f"product-boundary finite footprint differs: {key}")
    if settings.get("bias_correction") != "none" or report.get("automatic_bias_correction_applied") is not False:
        raise ValueError("this analysis preserves the ORAS5 fields without automatic bias correction")
    return report


def ensure_output_target(path: Path, *, overwrite: bool) -> None:
    resolved = Path(path).resolve()
    if resolved != PACKAGE_ROOT and PACKAGE_ROOT not in resolved.parents:
        raise ValueError(f"derived output must remain under {PACKAGE_ROOT}: {resolved}")
    if resolved.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite existing output: {resolved}")
    resolved.parent.mkdir(parents=True, exist_ok=True)


def atomic_replace(temp_path: Path, final_path: Path) -> None:
    temp_path = Path(temp_path).resolve()
    final_path = Path(final_path).resolve()
    if temp_path.parent != final_path.parent:
        raise ValueError("atomic replacement requires the same parent directory")
    temp_path.replace(final_path)


def write_json(path: Path, payload: Any, *, overwrite: bool = True) -> None:
    ensure_output_target(path, overwrite=overwrite)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    atomic_replace(temp, path)
