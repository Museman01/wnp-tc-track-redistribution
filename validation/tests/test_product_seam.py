"""Synthetic seam-diagnostic tests."""

from __future__ import annotations

from core.paths import PROJECT as _LAYOUT_ROOT

import sys
import unittest
from pathlib import Path

import numpy as np

SRC = _LAYOUT_ROOT / 'analysis/ocean'
sys.path.insert(0, str(SRC))
from analysis.ocean.ohc_product_boundary import empirical_absolute_percentile, robust_z  # noqa: E402
from analysis.ocean.ohc_common import coordinate_grids_equivalent  # noqa: E402


class TestProductSeam(unittest.TestCase):
    def test_reference_center_has_zero_score(self) -> None:
        self.assertAlmostEqual(robust_z(2.0, np.array([1.0, 2.0, 3.0])), 0.0)

    def test_large_jump_is_flaggable(self) -> None:
        score = robust_z(20.0, np.array([-1.0, 0.0, 1.0, 0.5, -0.5]))
        self.assertGreater(abs(score), 5.0)

    def test_constant_reference_is_explicit(self) -> None:
        self.assertTrue(np.isinf(robust_z(1.0, np.zeros(5))))

    def test_empirical_absolute_percentile_is_rank_based(self) -> None:
        percentile = empirical_absolute_percentile(2.5, np.array([-1.0, 2.0, 3.0]))
        self.assertAlmostEqual(percentile, 0.75)

    def test_coordinate_quantization_within_tolerance_is_equivalent(self) -> None:
        reference = np.array([[100.0, 100.25]], dtype="float32")
        current = reference.astype("float64") + np.array([[3.0e-5, -3.0e-5]])
        self.assertTrue(coordinate_grids_equivalent(current, reference, 5.0e-5))
        self.assertFalse(coordinate_grids_equivalent(current, reference, 2.0e-5))


if __name__ == "__main__":
    unittest.main(verbosity=2)
