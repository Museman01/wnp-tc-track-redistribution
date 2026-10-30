"""Known-pass and known-fail checks for product-boundary input validation."""

from core.paths import PROJECT as _LAYOUT_ROOT
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(_LAYOUT_ROOT / 'analysis/ocean'))
from analysis.ocean import ohc_common as common
@pytest.fixture
def seam(tmp_path, monkeypatch):
    monkeypatch.setattr(common, 'PROJECT_ROOT', tmp_path)
    report = {
        'structural_errors': [], 'left_year': 2014, 'right_year': 2015,
        'coordinate_tolerance_degrees': 0.00005,
        'seam_max_abs_latitude_difference_degrees': 0.0,
        'seam_max_abs_longitude_difference_degrees': 0.0,
        'seam_changed_finite_any_cells': 0, 'seam_changed_finite_all_cells': 0,
        'automatic_bias_correction_applied': False,
    }
    config = {
        'ohc': {'coordinate_equality_atol_degrees': 0.00005},
        'seam_audit': {'output_json': 'seam.json', 'left_year': 2014,
                       'right_year': 2015, 'bias_correction': 'none'},
    }
    path = tmp_path / 'seam.json'
    path.write_text(json.dumps(report), encoding='utf-8')
    return config, report, path


def test_accepts_complete_numerical_diagnostics(seam):
    config, report, _ = seam
    assert common.validate_product_seam(config) == report


@pytest.mark.parametrize('field,value', [
    ('structural_errors', ['coordinate mismatch']),
    ('seam_max_abs_latitude_difference_degrees', 0.001),
    ('seam_max_abs_longitude_difference_degrees', float('nan')),
    ('seam_changed_finite_any_cells', 1),
    ('seam_changed_finite_all_cells', 1),
    ('left_year', 2013),
    ('coordinate_tolerance_degrees', 0.1),
    ('automatic_bias_correction_applied', True),
])
def test_rejects_incompatible_diagnostics(seam, field, value):
    config, report, path = seam
    report[field] = value
    path.write_text(json.dumps(report), encoding='utf-8')
    with pytest.raises(ValueError):
        common.validate_product_seam(config)


def test_missing_diagnostics_are_explicit(seam):
    config, _, _ = seam
    config['seam_audit']['output_json'] = 'missing.json'
    with pytest.raises(FileNotFoundError):
        common.validate_product_seam(config)


def test_unconfigured_correction_is_rejected(seam):
    config, _, _ = seam
    config['seam_audit']['bias_correction'] = 'automatic'
    with pytest.raises(ValueError):
        common.validate_product_seam(config)
