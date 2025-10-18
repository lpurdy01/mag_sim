"""Generate three-phase stator scenarios with configurable resolution profiles."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple


DEFAULT_PROFILES: Dict[str, Dict[str, float | int]] = {
    "ci": {
        "nx": 65,
        "ny": 65,
        "frames_per_cycle": 12,
        "cycles": 1,
        "slot_depth": 0.010,
        "slot_width_deg": 18.0,
    },
    "hires": {
        "nx": 401,
        "ny": 401,
        "frames_per_cycle": 120,
        "cycles": 3,
        "slot_depth": 0.010,
        "slot_width_deg": 12.0,
    },
}

def build_circle(radius: float, segments: int) -> List[List[float]]:
    return [
        [radius * math.cos(2.0 * math.pi * i / segments), radius * math.sin(2.0 * math.pi * i / segments)]
        for i in range(segments)
    ]


def build_slot_polygon(
    center_angle_rad: float,
    inner_radius: float,
    outer_radius: float,
    width_deg: float,
) -> List[List[float]]:
    half_width = math.radians(width_deg) * 0.5
    angles = [center_angle_rad - half_width, center_angle_rad + half_width]
    points = [
        [inner_radius * math.cos(angles[0]), inner_radius * math.sin(angles[0])],
        [outer_radius * math.cos(angles[0]), outer_radius * math.sin(angles[0])],
        [outer_radius * math.cos(angles[1]), outer_radius * math.sin(angles[1])],
        [inner_radius * math.cos(angles[1]), inner_radius * math.sin(angles[1])],
    ]
    return points


def compute_phase_currents(current_amp: float, omega: float, time: float) -> Tuple[float, float, float]:
    theta = omega * time
    ia = current_amp * math.sin(theta)
    ib = current_amp * math.sin(theta - 2.0 * math.pi / 3.0)
    ic = current_amp * math.sin(theta + 2.0 * math.pi / 3.0)
    return ia, ib, ic


def generate_scenario(profile: str, output_path: Path) -> None:
    if profile not in DEFAULT_PROFILES:
        raise ValueError(f"Unknown profile '{profile}'. Available profiles: {', '.join(DEFAULT_PROFILES)}")

    profile_cfg = DEFAULT_PROFILES[profile]
    nx = int(profile_cfg["nx"])
    ny = int(profile_cfg["ny"])
    frames_per_cycle = int(profile_cfg["frames_per_cycle"])
    cycles = int(profile_cfg["cycles"])
    slot_depth = float(profile_cfg["slot_depth"])
    slot_width_deg = float(profile_cfg["slot_width_deg"])

    domain_size = 0.14
    r_in = 0.040
    r_out = 0.055
    bore_fraction = 0.6
    current_amp = 30.0
    slot_turns = 60.0
    slot_fill_fraction = 0.55
    electrical_hz = 60.0

    total_frames = frames_per_cycle * cycles
    dt = 1.0 / electrical_hz / frames_per_cycle
    omega = 2.0 * math.pi * electrical_hz

    outer_polygon = build_circle(r_out, 256)
    bore_polygon = build_circle(r_in, 192)

    slot_outer_radius = min(r_out - 0.002, r_in + slot_depth)
    slot_inner_radius = r_in

    slot_definitions = [
        ("A", +1.0, math.radians(0.0), "A_pos"),
        ("B", -1.0, math.radians(60.0), "B_neg"),
        ("C", +1.0, math.radians(120.0), "C_pos"),
        ("A", -1.0, math.radians(180.0), "A_neg"),
        ("B", +1.0, math.radians(240.0), "B_pos"),
        ("C", -1.0, math.radians(300.0), "C_neg"),
    ]

    slot_polygons: List[List[List[float]]] = []
    sources = []
    for phase, orientation, angle, slot_id in slot_definitions:
        polygon = build_slot_polygon(angle, slot_inner_radius, slot_outer_radius, slot_width_deg)
        slot_polygons.append(polygon)
        sources.append(
            {
                "type": "current_region",
                "id": slot_id,
                "label": f"Phase {phase}{'+' if orientation > 0 else '-'}",
                "phase": phase,
                "orientation": orientation,
                "vertices": polygon,
                "I": 0.0,
                "turns": slot_turns,
                "fill_fraction": slot_fill_fraction,
            }
        )

    regions: List[Dict[str, object]] = [
        {"type": "uniform", "material": "air"},
        {"type": "polygon", "material": "stator_steel", "vertices": outer_polygon},
    ]
    for polygon in slot_polygons:
        regions.append({"type": "polygon", "material": "air", "vertices": polygon})
    regions.append({"type": "polygon", "material": "air", "vertices": bore_polygon})

    bore_probe_polygon = build_circle(r_in * bore_fraction, 96)

    timeline = []
    for frame_idx in range(total_frames):
        time = frame_idx * dt
        ia, ib, ic = compute_phase_currents(current_amp, omega, time)
        timeline.append(
            {
                "t": time,
                "phase_currents": {
                    "A": ia,
                    "B": ib,
                    "C": ic,
                },
            }
        )

    scenario = {
        "version": "0.2",
        "units": "SI",
        "domain": {"Lx": domain_size, "Ly": domain_size, "nx": nx, "ny": ny},
        "boundary": {"type": "neumann"},
        "materials": [
            {"name": "air", "mu_r": 1.0},
            {"name": "stator_steel", "mu_r": 400.0},
        ],
        "regions": regions,
        "sources": sources,
        "timeline": timeline,
        "outputs": [
            {
                "type": "vtk_field_series",
                "id": "series",
                "dir": "outputs",
                "basename": "three_phase_frame",
                "quantities": ["B", "H"],
            },
            {
                "type": "polyline_outlines",
                "id": "outlines",
                "path": "outputs/three_phase_outlines.vtp",
            },
            {
                "type": "bore_avg_B",
                "id": "bore",
                "path": "outputs/bore_angle.csv",
                "vertices": bore_probe_polygon,
            },
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(scenario, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(DEFAULT_PROFILES.keys()), default="ci")
    parser.add_argument("--out", type=Path, required=True, help="Output scenario JSON path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_scenario(args.profile, args.out)


if __name__ == "__main__":
    main()
