#!/usr/bin/env python3
"""Smoke test for the rotor animation CLI."""

from __future__ import annotations

import copy
import csv
import json
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


def load_scenario(path: Path) -> tuple[dict, list[dict]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    timeline = data.get("timeline")
    if not timeline:
        raise SystemExit("Scenario timeline is missing for rotor animation smoke test")
    return data, timeline


def initial_rotor_angle_deg(spec: dict, rotor: str) -> float:
    mechanical = spec.get("mechanical", {}) or {}
    for entry in mechanical.get("rotors", []):
        if entry.get("name") == rotor:
            return float(entry.get("initial_angle_deg", 0.0))
    return 0.0


def synthesise_rotor_motion(times: list[float], initial_angle_deg: float) -> tuple[list[float], list[float]]:
    if not times:
        raise SystemExit("Timeline has no samples for rotor animation smoke test")
    angles_deg: list[float] = []
    omega_rad: list[float] = []
    for idx, current_time in enumerate(times):
        if idx == 0:
            angle = initial_angle_deg
            omega = 0.0
        else:
            dt = current_time - times[idx - 1]
            if dt <= 0.0:
                raise SystemExit("Scenario timeline is not strictly increasing")
            angle = initial_angle_deg + 6.0 * idx
            omega = math.radians(angle - angles_deg[-1]) / dt
        angles_deg.append(angle)
        omega_rad.append(omega)
    return angles_deg, omega_rad


def write_mechanical_trace(path: Path, rotor: str, times: list[float], angles_deg: list[float], omega_rad: list[float]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "rotor", "angle_deg", "omega_rad_s", "omega_rpm", "torque_Nm"])
        for t, angle, w in zip(times, angles_deg, omega_rad):
            rpm = w * 60.0 / (2.0 * math.pi)
            writer.writerow([t, rotor, angle, w, rpm, 0.05])


def write_circuit_trace(path: Path, times: list[float]) -> None:
    slot_ids = [
        "A_pos",
        "A_neg",
        "B_pos",
        "B_neg",
        "C_pos",
        "C_neg",
    ]
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


def write_timeline_with_rotor_angles(path: Path, spec: dict, rotor: str, angles_deg: list[float]) -> None:
    patched = copy.deepcopy(spec)
    timeline = patched.get("timeline", [])
    for frame, angle in zip(timeline, angles_deg):
        rotor_angles = frame.setdefault("rotor_angles", {})
        rotor_angles[rotor] = angle
    path.write_text(json.dumps(patched, indent=2), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    mechanical_csv = OUTPUT_DIR / "rotor_anim_mechanical.csv"
    circuit_csv = OUTPUT_DIR / "rotor_anim_circuit.csv"
    gif_path = OUTPUT_DIR / "rotor_anim.gif"
    png_path = OUTPUT_DIR / "rotor_anim.png"
    timeline_scenario_path = OUTPUT_DIR / "pm_motor_spinup_timeline.json"

    spec, timeline = load_scenario(SCENARIO)
    times = [float(frame.get("t", 0.0)) for frame in timeline]
    rotor_name = "pm_rotor"
    initial_angle = initial_rotor_angle_deg(spec, rotor_name)
    angles_deg, omega_rad = synthesise_rotor_motion(times, initial_angle)

    write_mechanical_trace(mechanical_csv, rotor_name, times, angles_deg, omega_rad)
    write_circuit_trace(circuit_csv, times)
    write_timeline_with_rotor_angles(timeline_scenario_path, spec, rotor_name, angles_deg)

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
    # the explicit traces. A patched copy of the scenario injects synthetic rotor
    # angles derived from the same motion used for the CSV traces so the helper
    # can operate without additional files.
    timeline_gif = OUTPUT_DIR / "rotor_anim_timeline.gif"
    timeline_png = OUTPUT_DIR / "rotor_anim_timeline.png"
    timeline_cmd = [
        sys.executable,
        str(ROOT / "python" / "generate_rotor_animation.py"),
        "--scenario",
        str(timeline_scenario_path),
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
