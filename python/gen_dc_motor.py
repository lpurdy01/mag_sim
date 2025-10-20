"""Generate DC motor scenarios with commutated armature support."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List

Number = float

DEFAULT_PROFILES: Dict[str, Dict[str, Number]] = {
    "ci": {
        "nx": 65,
        "ny": 65,
        "frames_per_cycle": 12,
        "cycles": 2,
        "frame_dt": 1.0 / 720.0,
        "outer_segments": 48,
        "bore_segments": 36,
        "rotor_segments": 24,
        "torque_segments": 40,
        "quantize_digits": 3,
    },
    "hires": {
        "nx": 401,
        "ny": 401,
        "frames_per_cycle": 120,
        "cycles": 6,
        "frame_dt": 1.0 / 2400.0,
        "outer_segments": 360,
        "bore_segments": 288,
        "rotor_segments": 216,
        "torque_segments": 288,
        "quantize_digits": 6,
    },
}

VALID_MODES = {"locked", "spinup", "commutator_test"}


def build_circle(radius: float, segments: int) -> List[List[float]]:
    return [
        [radius * math.cos(2.0 * math.pi * i / segments), radius * math.sin(2.0 * math.pi * i / segments)]
        for i in range(segments)
    ]


def build_arc_segment(center_deg: float, inner_radius: float, outer_radius: float, width_deg: float) -> List[List[float]]:
    half_width = math.radians(width_deg * 0.5)
    center = math.radians(center_deg)
    start = center - half_width
    end = center + half_width
    return [
        [inner_radius * math.cos(start), inner_radius * math.sin(start)],
        [outer_radius * math.cos(start), outer_radius * math.sin(start)],
        [outer_radius * math.cos(end), outer_radius * math.sin(end)],
        [inner_radius * math.cos(end), inner_radius * math.sin(end)],
    ]


def quantize_numbers(value, ndigits: int):
    if isinstance(value, float):
        quantized = round(value, ndigits)
        if quantized == -0.0:
            quantized = 0.0
        return quantized
    if isinstance(value, list):
        return [quantize_numbers(item, ndigits) for item in value]
    if isinstance(value, dict):
        return {key: quantize_numbers(val, ndigits) for key, val in value.items()}
    return value


def make_current_region(
    *,
    region_id: str,
    label: str,
    phase: str,
    orientation: float,
    turns: float,
    fill_fraction: float,
    polygon: List[List[float]],
) -> Dict[str, object]:
    return {
        "type": "current_region",
        "id": region_id,
        "label": label,
        "phase": phase,
        "orientation": orientation,
        "vertices": polygon,
        "I": 0.0,
        "turns": turns,
        "fill_fraction": fill_fraction,
    }


def build_base_scenario(profile: str, mode: str) -> Dict[str, object]:
    if profile not in DEFAULT_PROFILES:
        raise ValueError(f"Unknown profile '{profile}'. Available profiles: {', '.join(DEFAULT_PROFILES)}")
    if mode not in VALID_MODES:
        raise ValueError(f"Unknown mode '{mode}'. Valid options: {', '.join(sorted(VALID_MODES))}")

    cfg = DEFAULT_PROFILES[profile]
    nx = int(cfg["nx"])
    ny = int(cfg["ny"])
    frame_dt = float(cfg["frame_dt"])
    frames_per_cycle = int(cfg["frames_per_cycle"])
    cycles = int(cfg["cycles"])
    total_frames = max(1, frames_per_cycle * cycles)

    outer_segments = int(cfg["outer_segments"])
    bore_segments = int(cfg["bore_segments"])
    rotor_segments = int(cfg["rotor_segments"])
    torque_segments = int(cfg["torque_segments"])

    domain_radius = 0.11
    domain_size = domain_radius * 2.0

    stator_outer_radius = 0.055
    bore_radius = 0.020
    rotor_radius = 0.015
    rotor_slot_inner = 0.009
    rotor_slot_outer = 0.014
    field_slot_inner = bore_radius + 0.0015
    field_slot_outer = bore_radius + 0.008

    outer_polygon = build_circle(stator_outer_radius, outer_segments)
    bore_polygon = build_circle(bore_radius, bore_segments)
    rotor_polygon = build_circle(rotor_radius, rotor_segments)
    torque_loop = build_circle(bore_radius + 0.002, torque_segments)

    field_turns = 220.0
    field_fill = 0.58
    field_width_deg = 70.0
    north_polygon = build_arc_segment(90.0, field_slot_inner, field_slot_outer, field_width_deg)
    south_polygon = build_arc_segment(270.0, field_slot_inner, field_slot_outer, field_width_deg)

    arm_turns = 110.0
    arm_fill = 0.52
    arm_width_deg = 80.0
    arm_polygon_a = build_arc_segment(0.0, rotor_slot_inner, rotor_slot_outer, arm_width_deg)
    arm_polygon_b = build_arc_segment(180.0, rotor_slot_inner, rotor_slot_outer, arm_width_deg)

    regions: List[Dict[str, object]] = [{"type": "uniform", "material": "air"}]
    regions.append({"type": "polygon", "material": "stator_steel", "vertices": outer_polygon})
    regions.append({"type": "polygon", "material": "air", "vertices": north_polygon})
    regions.append({"type": "polygon", "material": "air", "vertices": south_polygon})
    regions.append({"type": "polygon", "material": "air", "vertices": bore_polygon})
    rotor_polygon_index = len(regions)
    regions.append({"type": "polygon", "material": "rotor_steel", "vertices": rotor_polygon})
    regions.append({"type": "polygon", "material": "air", "vertices": arm_polygon_a})
    regions.append({"type": "polygon", "material": "air", "vertices": arm_polygon_b})

    sources = [
        make_current_region(
            region_id="field_north",
            label="Field North",
            phase="field",
            orientation=1.0,
            turns=field_turns,
            fill_fraction=field_fill,
            polygon=north_polygon,
        ),
        make_current_region(
            region_id="field_south",
            label="Field South",
            phase="field",
            orientation=-1.0,
            turns=field_turns,
            fill_fraction=field_fill,
            polygon=south_polygon,
        ),
        make_current_region(
            region_id="armature_a",
            label="Armature A",
            phase="armature",
            orientation=1.0,
            turns=arm_turns,
            fill_fraction=arm_fill,
            polygon=arm_polygon_a,
        ),
        make_current_region(
            region_id="armature_b",
            label="Armature B",
            phase="armature",
            orientation=-1.0,
            turns=arm_turns,
            fill_fraction=arm_fill,
            polygon=arm_polygon_b,
        ),
    ]

    rotors = [
        {
            "name": "dc_rotor",
            "pivot": [0.0, 0.0],
            "polygons": [rotor_polygon_index],
            "initial_angle_deg": -25.0,
        }
    ]

    circuits = [
        {
            "id": "dc_drive",
            "nodes": ["gnd", "field_pos", "arm_pos"],
            "elements": [
                {"type": "voltage_source", "id": "field_V", "nodes": ["field_pos", "gnd"], "value": 18.0},
                {"type": "resistor", "id": "field_R", "nodes": ["field_pos", "gnd"], "resistance": 1.6},
                {"type": "inductor", "id": "field_L", "nodes": ["field_pos", "gnd"], "inductance": 0.08, "initial_current": 0.0},
                {"type": "voltage_source", "id": "arm_V", "nodes": ["arm_pos", "gnd"], "value": 12.0},
                {"type": "resistor", "id": "arm_R", "nodes": ["arm_pos", "gnd"], "resistance": 0.75},
                {"type": "inductor", "id": "arm_L", "nodes": ["arm_pos", "gnd"], "inductance": 0.032, "initial_current": 0.0},
                {"type": "coil_link", "id": "field_coil_north", "inductor": "field_L", "region": "field_north", "turns": field_turns},
                {"type": "coil_link", "id": "field_coil_south", "inductor": "field_L", "region": "field_south", "turns": field_turns},
                {
                    "type": "coil_link",
                    "id": "arm_coil_a",
                    "inductor": "arm_L",
                    "region": "armature_a",
                    "turns": arm_turns,
                    "commutator": {
                        "rotor": "dc_rotor",
                        "default_orientation": 1.0,
                        "segments": [
                            {"start_deg": -90.0, "end_deg": 90.0, "orientation": 1.0},
                            {"start_deg": 90.0, "end_deg": 270.0, "orientation": -1.0},
                        ],
                    },
                },
                {
                    "type": "coil_link",
                    "id": "arm_coil_b",
                    "inductor": "arm_L",
                    "region": "armature_b",
                    "turns": arm_turns,
                    "commutator": {
                        "rotor": "dc_rotor",
                        "default_orientation": 1.0,
                        "segments": [
                            {"start_deg": -90.0, "end_deg": 90.0, "orientation": 1.0},
                            {"start_deg": 90.0, "end_deg": 270.0, "orientation": -1.0},
                        ],
                    },
                },
            ],
        }
    ]

    scenario: Dict[str, object] = {
        "version": "0.2",
        "units": "SI",
        "domain": {"Lx": domain_size, "Ly": domain_size, "nx": nx, "ny": ny},
        "boundary": {"type": "neumann"},
        "materials": [
            {"name": "air", "mu_r": 1.0},
            {"name": "stator_steel", "mu_r": 1800.0},
            {"name": "rotor_steel", "mu_r": 1200.0},
        ],
        "regions": regions,
        "rotors": rotors,
        "sources": sources,
        "circuits": circuits,
    }

    outputs = [
        {"type": "vtk_field_series", "id": "dc_motor_series", "dir": "outputs", "basename": "dc_motor_frame", "quantities": ["B", "H"]},
        {"type": "polyline_outlines", "id": "dc_motor_outlines", "path": "outputs/dc_motor_outlines.vtp"},
        {
            "type": "circuit_trace",
            "id": "dc_motor_currents",
            "path": "outputs/dc_motor_currents.csv",
            "circuits": ["dc_drive"],
        },
        {
            "type": "probe",
            "id": "dc_torque",
            "probe_type": "torque",
            "method": "stress_tensor",
            "loop": {"type": "polygon", "vertices": torque_loop},
            "path": "outputs/dc_motor_torque.csv",
        },
    ]

    scenario["outputs"] = outputs

    if mode == "spinup":
        scenario["mechanical"] = {
            "rotors": [
                {
                    "name": "dc_rotor",
                    "inertia": 0.00065,
                    "damping": 0.000045,
                    "load_torque": 0.08,
                    "initial_angle_deg": -25.0,
                    "initial_speed_rpm": 0.0,
                    "torque_probe": "dc_torque",
                }
            ]
        }
        scenario.setdefault("outputs", []).append(
            {
                "type": "mechanical_trace",
                "id": "dc_motor_mechanical",
                "path": "outputs/dc_motor_mechanical.csv",
                "rotors": ["dc_rotor"],
            }
        )
        scenario["timeline"] = [{"t": frame_idx * frame_dt} for frame_idx in range(total_frames)]
    elif mode == "locked":
        scenario["mechanical"] = {
            "rotors": [
                {
                    "name": "dc_rotor",
                    "inertia": 0.00065,
                    "damping": 0.00005,
                    "load_torque": 0.1,
                    "initial_angle_deg": -20.0,
                    "torque_probe": "dc_torque",
                }
            ]
        }
        scenario["timeline"] = [{"t": frame_idx * frame_dt} for frame_idx in range(total_frames)]
    else:  # commutator_test
        angles = [-135.0, -45.0, 45.0, 135.0]
        scenario["timeline"] = [
            {"t": idx * frame_dt, "rotor_angles": {"dc_rotor": angle}} for idx, angle in enumerate(angles)
        ]
        elements = []
        for element in circuits[0]["elements"]:
            if element.get("id") == "arm_L":
                updated = dict(element)
                updated["initial_current"] = 8.0
                elements.append(updated)
            else:
                elements.append(element)
        circuits[0]["elements"] = elements

    return scenario


def generate_scenario(profile: str, mode: str, output_path: Path) -> None:
    scenario = build_base_scenario(profile, mode)
    quantize_digits = int(DEFAULT_PROFILES[profile].get("quantize_digits", 6))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(quantize_numbers(scenario, quantize_digits), indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="ci", choices=DEFAULT_PROFILES.keys())
    parser.add_argument("--mode", default="spinup", choices=sorted(VALID_MODES))
    parser.add_argument("--out", required=True, type=Path, help="Output scenario path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_scenario(args.profile, args.mode, args.out)


if __name__ == "__main__":
    main()
