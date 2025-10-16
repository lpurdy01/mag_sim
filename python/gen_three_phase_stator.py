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
        "load_angle_deg": 15.0,
    },
    "hires": {
        "nx": 401,
        "ny": 401,
        "frames_per_cycle": 120,
        "cycles": 3,
        "slot_depth": 0.010,
        "slot_width_deg": 12.0,
        "load_angle_deg": 10.0,
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
    load_angle_deg = float(profile_cfg.get("load_angle_deg", 15.0))

    domain_size = 0.14
    r_in = 0.040
    r_out = 0.055
    bore_fraction = 0.6
    current_amp = 200.0
    electrical_hz = 60.0
    magnet_strength = 4.0e5
    phase_resistance = 0.4
    phase_inductance = 0.012
    coil_turns = 120.0
    line_voltage_rms = 325.0  # peak phase voltage ~230 RMS line-to-neutral
    rotor_inertia = 0.0025
    rotor_damping = 0.0002
    load_torque = 8.0

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
    phase_slot_map: Dict[str, Dict[str, List[List[float]]]] = {"A": {}, "B": {}, "C": {}}
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
            }
        )
        phase_slot_map.setdefault(phase, {})["pos" if orientation > 0 else "neg"] = polygon

    regions: List[Dict[str, object]] = [
        {"type": "uniform", "material": "air"},
        {"type": "polygon", "material": "stator_steel", "vertices": outer_polygon},
    ]
    next_polygon_index = 1  # stator polygon inserted above
    for polygon in slot_polygons:
        regions.append({"type": "polygon", "material": "air", "vertices": polygon})
        next_polygon_index += 1
    regions.append({"type": "polygon", "material": "air", "vertices": bore_polygon})
    next_polygon_index += 1

    rotor_radius = r_in * 0.55
    rotor_polygon = build_circle(rotor_radius, 160)
    rotor_polygon_index = next_polygon_index
    regions.append({"type": "polygon", "material": "rotor_steel", "vertices": rotor_polygon})
    next_polygon_index += 1

    magnet_half_width = rotor_radius * 0.85
    magnet_height = rotor_radius * 0.55
    magnet_polygon = [
        [-magnet_half_width, -magnet_height],
        [magnet_half_width, -magnet_height],
        [magnet_half_width, magnet_height],
        [-magnet_half_width, magnet_height],
    ]

    magnet_regions = [
        {
            "type": "polygon",
            "vertices": magnet_polygon,
            "magnetization": [magnet_strength, 0.0],
        }
    ]

    rotors = [
        {
            "name": "pm_rotor",
            "pivot": [0.0, 0.0],
            "polygons": [rotor_polygon_index],
            "magnets": [0],
        }
    ]

    bore_probe_polygon = build_circle(r_in * bore_fraction, 96)
    torque_loop_radius = rotor_radius + 0.004
    torque_probe_polygon = build_circle(torque_loop_radius, 128)

    timeline = []
    phase_voltage_amp = line_voltage_rms
    for frame_idx in range(total_frames):
        time = frame_idx * dt
        ia, ib, ic = compute_phase_currents(current_amp, omega, time)
        va = phase_voltage_amp * math.sin(omega * time)
        vb = phase_voltage_amp * math.sin(omega * time - 2.0 * math.pi / 3.0)
        vc = phase_voltage_amp * math.sin(omega * time + 2.0 * math.pi / 3.0)
        electrical_angle_deg = math.degrees(omega * time)
        rotor_angle_deg = electrical_angle_deg - load_angle_deg
        timeline.append(
            {
                "t": time,
                "phase_currents": {
                    "A": ia,
                    "B": ib,
                    "C": ic,
                },
                "voltage_sources": [
                    {"circuit": "stator_three_phase", "source": "Va", "value": va},
                    {"circuit": "stator_three_phase", "source": "Vb", "value": vb},
                    {"circuit": "stator_three_phase", "source": "Vc", "value": vc},
                ],
                "rotor_angles": {"pm_rotor": rotor_angle_deg},
            }
        )

    circuits = [
        {
            "id": "stator_three_phase",
            "nodes": ["neutral", "phase_a", "phase_b", "phase_c"],
            "elements": [
                {"type": "resistor", "id": "R_a", "nodes": ["phase_a", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_a", "nodes": ["phase_a", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Va", "nodes": ["phase_a", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_a_pos", "inductor": "L_a", "region": "A_pos", "turns": coil_turns},
                {"type": "coil_link", "id": "coil_a_neg", "inductor": "L_a", "region": "A_neg", "turns": coil_turns},
                {"type": "resistor", "id": "R_b", "nodes": ["phase_b", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_b", "nodes": ["phase_b", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Vb", "nodes": ["phase_b", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_b_pos", "inductor": "L_b", "region": "B_pos", "turns": coil_turns},
                {"type": "coil_link", "id": "coil_b_neg", "inductor": "L_b", "region": "B_neg", "turns": coil_turns},
                {"type": "resistor", "id": "R_c", "nodes": ["phase_c", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_c", "nodes": ["phase_c", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Vc", "nodes": ["phase_c", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_c_pos", "inductor": "L_c", "region": "C_pos", "turns": coil_turns},
                {"type": "coil_link", "id": "coil_c_neg", "inductor": "L_c", "region": "C_neg", "turns": coil_turns},
            ],
        }
    ]

    scenario = {
        "version": "0.2",
        "units": "SI",
        "domain": {"Lx": domain_size, "Ly": domain_size, "nx": nx, "ny": ny},
        "boundary": {"type": "neumann"},
        "materials": [
            {"name": "air", "mu_r": 1.0},
            {"name": "stator_steel", "mu_r": 400.0},
            {"name": "rotor_steel", "mu_r": 800.0},
        ],
        "regions": regions,
        "magnet_regions": magnet_regions,
        "rotors": rotors,
        "sources": sources,
        "circuits": circuits,
        "mechanical": {
            "rotors": [
                {
                    "name": "pm_rotor",
                    "inertia": rotor_inertia,
                    "damping": rotor_damping,
                    "load_torque": load_torque,
                    "initial_angle_deg": -load_angle_deg,
                    "initial_speed_rpm": 0.0,
                    "torque_probe": "rotor_torque",
                }
            ]
        },
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
            {
                "type": "probe",
                "id": "rotor_torque",
                "probe_type": "torque",
                "method": "stress_tensor",
                "loop": {"type": "polygon", "vertices": torque_probe_polygon},
                "path": "outputs/rotor_torque.csv",
            },
        ],
    }

    back_emf_outputs = []
    for phase, slots in phase_slot_map.items():
        pos_polygon = slots.get("pos")
        if not pos_polygon:
            continue
        back_emf_outputs.append(
            {
                "type": "back_emf_probe",
                "id": f"phase_{phase.lower()}_emf",
                "component": "Bmag",
                "region": {"type": "polygon", "vertices": pos_polygon},
                "path": f"outputs/phase_{phase.lower()}_emf.csv",
            }
        )

    scenario["outputs"].extend(back_emf_outputs)

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
