"""Sanity checks for the three-phase stator rotating field simulation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from animate_three_phase import (  # reuse data loaders
    compute_bore_series,
    extract_bore_polygon,
    load_pvd,
    load_scenario,
)


def compute_r_squared(x: np.ndarray, y: np.ndarray) -> float:
    coeffs = np.polyfit(x, y, 1)
    fit = np.polyval(coeffs, x)
    ss_res = np.sum((y - fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def check_rotating_field(pvd_path: str, scenario_path: str) -> int:
    scenario = load_scenario(Path(scenario_path))
    datasets = load_pvd(Path(pvd_path))
    if not datasets:
        print("No datasets found in PVD file")
        return 1
    bore_polygon = extract_bore_polygon(scenario)
    times, bx, by = compute_bore_series(datasets, bore_polygon)
    magnitude = np.hypot(bx, by)
    angles = np.unwrap(np.arctan2(by, bx))
    diffs = np.diff(angles)
    direction = 1.0 if np.mean(diffs) >= 0.0 else -1.0
    projected = direction * diffs
    min_step = float(np.min(projected)) if projected.size else 0.0
    if min_step <= -2e-3:
        print(
            "Rotating field check FAILED: angular steps change sign (min projected Δ={:.4e})".format(
                min_step
            )
        )
        return 1
    r_squared = compute_r_squared(times, angles)
    if r_squared < 0.95:
        print(f"Rotating field check FAILED: R^2={r_squared:.3f} < 0.95")
        return 1
    if np.min(magnitude) < 1e-4:
        print("Rotating field check FAILED: |B| collapsed below threshold")
        return 1
    print(
        "Rotating field check PASSED: angle monotonic, R^2={:.3f}, |B|min={:.4f}, min step={:.4e}".format(
            r_squared, np.min(magnitude), min_step
        )
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pvd", required=True, help="VTK time-series PVD file")
    parser.add_argument("--scenario", required=True, help="Scenario JSON path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    exit(check_rotating_field(args.pvd, args.scenario))


if __name__ == "__main__":
    main()
