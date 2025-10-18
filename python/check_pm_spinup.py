#!/usr/bin/env python3
"""Validate mechanical spin-up traces emitted by synchronous or induction demos."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import List, Tuple


def load_mechanical_trace(path: Path, rotor: str) -> Tuple[List[float], List[float], List[float]]:
    times: List[float] = []
    angles: List[float] = []
    speeds: List[float] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("rotor") != rotor:
                continue
            times.append(float(row["time_s"]))
            angles.append(float(row["angle_deg"]))
            speeds.append(float(row["omega_rad_s"]))
    return times, angles, speeds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mechanical", required=True, type=Path, help="Path to the mechanical trace CSV")
    parser.add_argument("--scenario", required=True, type=Path, help="Scenario JSON used for the run")
    parser.add_argument("--rotor", default="pm_rotor", help="Rotor name to analyse")
    parser.add_argument("--min-angle-rise", type=float, default=15.0, help="Minimum required angle increase (deg)")
    parser.add_argument("--min-speed-rise", type=float, default=10.0, help="Minimum required speed rise (rad/s)")
    parser.add_argument("--max-backstep", type=float, default=1.0, help="Allowed negative angle step (deg)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.mechanical.exists():
        raise SystemExit(f"Mechanical trace '{args.mechanical}' not found")
    if not args.scenario.exists():
        raise SystemExit(f"Scenario '{args.scenario}' not found")

    scenario = json.loads(args.scenario.read_text(encoding="utf-8"))
    rotors = [rotor.get("name", "") for rotor in scenario.get("mechanical", {}).get("rotors", [])]
    if args.rotor not in rotors:
        raise SystemExit(f"Rotor '{args.rotor}' not declared in scenario mechanical section; available: {rotors}")

    times, angles, speeds = load_mechanical_trace(args.mechanical, args.rotor)
    if len(times) < 3:
        raise SystemExit(f"Mechanical trace for rotor '{args.rotor}' contains too few samples ({len(times)})")

    paired = sorted(zip(times, angles, speeds), key=lambda entry: entry[0])
    times = [p[0] for p in paired]
    angles = [p[1] for p in paired]
    speeds = [p[2] for p in paired]

    angle_delta = angles[-1] - angles[0]
    direction = 1.0 if angle_delta >= 0.0 else -1.0
    progress = [direction * (angles[i + 1] - angles[i]) for i in range(len(angles) - 1)]
    min_progress = min(progress)
    if min_progress < -abs(args.max_backstep):
        raise SystemExit(
            "Rotor angle regressed by "
            f"{(-direction * min_progress):.3f} deg which exceeds allowed backstep {abs(args.max_backstep):.3f} deg"
        )

    angle_magnitude = abs(angle_delta)
    if angle_magnitude < args.min_angle_rise:
        raise SystemExit(
            f"Insufficient angle change for rotor '{args.rotor}': |Δθ|={angle_magnitude:.2f} deg < {args.min_angle_rise:.2f} deg"
        )

    speed_delta = speeds[-1] - speeds[0]
    speed_magnitude = abs(speeds[-1]) - abs(speeds[0])
    if speed_magnitude < args.min_speed_rise:
        raise SystemExit(
            f"Insufficient speed change for rotor '{args.rotor}': |Δω|={speed_magnitude:.2f} rad/s < {args.min_speed_rise:.2f} rad/s"
        )
    if abs(speeds[-1]) <= 0.0:
        raise SystemExit(f"Final rotor speed magnitude is zero ({speeds[-1]:.3f} rad/s)")

    rpm = abs(speeds[-1]) * 60.0 / (2.0 * math.pi)
    direction_label = "increased" if direction >= 0.0 else "decreased"
    print(
        f"Spin-up check passed: Δθ={angle_delta:.2f} deg (|Δθ|={angle_magnitude:.2f} deg, {direction_label}), "
        f"Δω={speed_delta:.2f} rad/s (|Δω|={speed_magnitude:.2f} rad/s), final |ω|={abs(speeds[-1]):.2f} rad/s "
        f"({rpm:.1f} rpm) over {len(times)} samples"
    )


if __name__ == "__main__":
    main()
