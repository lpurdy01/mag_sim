#!/usr/bin/env python3
"""Smoke test for the rotor animation CLI."""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from pathlib import Path

try:  # pragma: no cover - optional dependency guard
    import numpy  # noqa: F401
    import matplotlib  # noqa: F401
except ModuleNotFoundError as exc:  # pragma: no cover - graceful skip
    print(f"Skipping rotor animation smoke test: missing dependency ({exc})")
    sys.exit(0)


ROOT = Path(__file__).resolve().parent.parent
SCENARIO = ROOT / "inputs" / "tests" / "pm_motor_spinup_test.json"
OUTPUT_DIR = ROOT / "outputs"


def write_mechanical_trace(path: Path, rotor: str) -> None:
    times = [0.0, 0.001, 0.002, 0.003, 0.004]
    angles_deg = [0.0, 2.0, 4.5, 7.0, 9.5]
    omega = [math.radians(angles_deg[i] - angles_deg[i - 1]) / (times[i] - times[i - 1]) if i > 0 else 0.0 for i in range(len(times))]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "rotor", "angle_deg", "omega_rad_s", "omega_rpm", "torque_Nm"])
        for t, angle, w in zip(times, angles_deg, omega):
            rpm = w * 60.0 / (2.0 * math.pi)
            writer.writerow([t, rotor, angle, w, rpm, 0.05])


def write_circuit_trace(path: Path) -> None:
    slot_ids = [
        "A_pos",
        "A_neg",
        "B_pos",
        "B_neg",
        "C_pos",
        "C_neg",
    ]
    times = [0.0, 0.001, 0.002, 0.003, 0.004]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_s",
                "circuit",
                "coil",
                "region",
                "turns",
                "orientation",
                "current_A",
                "effective_current_A",
                "ampere_turns",
                "back_emf_V",
            ]
        )
        for t in times:
            phase_angle = 2.0 * math.pi * 60.0 * t
            for idx, slot in enumerate(slot_ids):
                phase_offset = idx * (math.pi / 3.0)
                current = 10.0 * math.sin(phase_angle - phase_offset)
                ampere_turns = current * 120.0
                writer.writerow([t, "stator_three_phase", slot, slot, 120.0, 1.0, current, current, ampere_turns, 0.0])


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mechanical_csv = OUTPUT_DIR / "rotor_anim_mechanical.csv"
    circuit_csv = OUTPUT_DIR / "rotor_anim_circuit.csv"
    gif_path = OUTPUT_DIR / "rotor_anim.gif"
    png_path = OUTPUT_DIR / "rotor_anim.png"

    write_mechanical_trace(mechanical_csv, "pm_rotor")
    write_circuit_trace(circuit_csv)

    cmd = [
        sys.executable,
        str(ROOT / "python" / "generate_rotor_animation.py"),
        "--scenario",
        str(SCENARIO),
        "--rotor",
        "pm_rotor",
        "--mechanical",
        str(mechanical_csv),
        "--mechanical-rotor",
        "pm_rotor",
        "--circuit-trace",
        str(circuit_csv),
        "--gif",
        str(gif_path),
        "--frame-png",
        str(png_path),
        "--max-frames",
        "5",
        "--fps",
        "10",
    ]
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    subprocess.run(cmd, check=True, cwd=ROOT, env=env)

    for artifact in (gif_path, png_path):
        if not artifact.exists() or artifact.stat().st_size == 0:
            raise SystemExit(f"Expected artifact {artifact} was not created")

    # Exercise the timeline fallbacks (rotor angles + phase currents) by omitting
    # the explicit traces. The PM spin-up scenario ships both series, so the
    # helper should succeed without additional CSV inputs.
    timeline_gif = OUTPUT_DIR / "rotor_anim_timeline.gif"
    timeline_png = OUTPUT_DIR / "rotor_anim_timeline.png"
    timeline_cmd = [
        sys.executable,
        str(ROOT / "python" / "generate_rotor_animation.py"),
        "--scenario",
        str(SCENARIO),
        "--rotor",
        "pm_rotor",
        "--gif",
        str(timeline_gif),
        "--frame-png",
        str(timeline_png),
        "--max-frames",
        "6",
        "--fps",
        "12",
    ]
    subprocess.run(timeline_cmd, check=True, cwd=ROOT, env=env)

    for artifact in (timeline_gif, timeline_png):
        if not artifact.exists() or artifact.stat().st_size == 0:
            raise SystemExit(f"Expected artifact {artifact} was not created")

    print("Rotor animation smoke test produced:")
    print(f"  GIF: {gif_path}")
    print(f"  PNG: {png_path}")
    print(f"  Timeline GIF: {timeline_gif}")
    print(f"  Timeline PNG: {timeline_png}")


if __name__ == "__main__":
    main()
