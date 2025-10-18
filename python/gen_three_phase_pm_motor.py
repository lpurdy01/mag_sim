"""Generate three-phase permanent-magnet motor scenarios with co-simulation data."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


Number = Union[int, float]


DEFAULT_PROFILES: Dict[str, Dict[str, Number]] = {
    "ci": {
        "nx": 65,
        "ny": 65,
        "frames_per_cycle": 12,
        "cycles": 1,
        "slot_depth": 0.010,
        "slot_width_deg": 18.0,
        "load_angle_deg": 15.0,
        "spinup_cycles": 1,
        "spinup_load_angle_deg": 30.0,
        "spinup_initial_angle_deg": -30.0,
        "spinup_load_torque": 0.12,
        "spinup_damping": 0.00003,
        "spinup_frames_per_cycle": 10,
        "outer_segments": 18,
        "bore_segments": 12,
        "rotor_segments": 10,
        "torque_segments": 24,
        "bore_probe_segments": 18,
        "quantize_digits": 3,
    },
    "hires": {
        "nx": 401,
        "ny": 401,
        "frames_per_cycle": 120,
        "cycles": 3,
        "slot_depth": 0.010,
        "slot_width_deg": 12.0,
        "load_angle_deg": 10.0,
        "spinup_cycles": 5,
        "spinup_load_angle_deg": 24.0,
        "spinup_initial_angle_deg": -24.0,
        "spinup_load_torque": 1.2,
        "spinup_damping": 0.0002,
        "spinup_frames_per_cycle": 100,
        "outer_segments": 384,
        "bore_segments": 288,
        "rotor_segments": 256,
        "torque_segments": 384,
        "bore_probe_segments": 288,
        "quantize_digits": 6,
    },
}

VALID_MODES = ("locked", "spinup")


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


def quantize_numbers(value: object, ndigits: int = 6) -> object:
    """Round floats recursively so CI fixtures stay compact."""

    if isinstance(value, float):
        quantized = round(value, ndigits)
        # Avoid signed zero noise in the JSON dumps.
        if quantized == 0.0:
            quantized = 0.0
        return quantized
    if isinstance(value, list):
        return [quantize_numbers(item, ndigits) for item in value]
    if isinstance(value, tuple):
        return [quantize_numbers(item, ndigits) for item in value]
    if isinstance(value, dict):
        return {key: quantize_numbers(val, ndigits) for key, val in value.items()}
    return value


def generate_scenario(
    profile: str,
    output_path: Path,
    *,
    mode: str = "locked",
    cycles_override: Optional[int] = None,
    frames_per_cycle_override: Optional[int] = None,
) -> None:
    if profile not in DEFAULT_PROFILES:
        raise ValueError(f"Unknown profile '{profile}'. Available profiles: {', '.join(DEFAULT_PROFILES)}")
    if mode not in VALID_MODES:
        raise ValueError(f"Unknown mode '{mode}'. Valid options: {', '.join(VALID_MODES)}")

    profile_cfg = DEFAULT_PROFILES[profile]
    nx = int(profile_cfg["nx"])
    ny = int(profile_cfg["ny"])

    frames_per_cycle = frames_per_cycle_override if frames_per_cycle_override is not None else int(
        profile_cfg["frames_per_cycle"]
    )
    base_cycles = int(profile_cfg["cycles"])
    cycles = cycles_override if cycles_override is not None else base_cycles

    slot_depth = float(profile_cfg["slot_depth"])
    slot_width_deg = float(profile_cfg["slot_width_deg"])
    base_load_angle_deg = float(profile_cfg.get("load_angle_deg", 15.0))
    load_angle_deg = base_load_angle_deg

    if mode == "spinup":
        if cycles_override is None:
            cycles = int(profile_cfg.get("spinup_cycles", base_cycles))
        load_angle_deg = float(profile_cfg.get("spinup_load_angle_deg", base_load_angle_deg))
        spinup_fpc = int(profile_cfg.get("spinup_frames_per_cycle", frames_per_cycle))
        if frames_per_cycle_override is None:
            frames_per_cycle = spinup_fpc

    domain_size = 0.14
    r_in = 0.040
    r_out = 0.055
    bore_fraction = 0.6
    current_amp = 35.0
    electrical_hz = 60.0
    magnet_strength = 1.0e5
    magnet_mu_r = 1.05
    phase_resistance = 0.4
    phase_inductance = 0.012
    phase_turns = 120.0
    turns_per_slot = phase_turns / 2.0
    slot_fill_fraction = 0.55
    line_voltage_rms = 20.0  # peak phase voltage driving the RL network
    rotor_inertia = 0.0008
    rotor_damping = 0.00005
    load_torque = 0.12

    if mode == "spinup":
        rotor_damping = float(profile_cfg.get("spinup_damping", rotor_damping))
        load_torque = float(profile_cfg.get("spinup_load_torque", load_torque))

    total_frames = frames_per_cycle * cycles
    dt = 1.0 / electrical_hz / frames_per_cycle
    omega = 2.0 * math.pi * electrical_hz

    outer_segments = int(profile_cfg.get("outer_segments", 256))
    bore_segments = int(profile_cfg.get("bore_segments", 192))
    rotor_segments = int(profile_cfg.get("rotor_segments", 160))
    torque_segments = int(profile_cfg.get("torque_segments", max(64, rotor_segments)))
    bore_probe_segments = int(profile_cfg.get("bore_probe_segments", max(64, bore_segments)))
    quantize_digits = int(profile_cfg.get("quantize_digits", 6))

    outer_polygon = build_circle(r_out, outer_segments)
    bore_polygon = build_circle(r_in, bore_segments)

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
                "turns": turns_per_slot,
                "fill_fraction": slot_fill_fraction,
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
    rotor_polygon = build_circle(rotor_radius, rotor_segments)
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

    regions.append({"type": "polygon", "material": "pm_magnet", "vertices": magnet_polygon})

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

    bore_probe_polygon = build_circle(r_in * bore_fraction, bore_probe_segments)
    torque_loop_radius = rotor_radius + 0.004
    torque_probe_polygon = build_circle(torque_loop_radius, torque_segments)

    prefix = "pm_motor_spinup" if mode == "spinup" else "pm_motor"
    torque_id = f"{prefix}_torque"

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
        frame_entry = {
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
        }
        if mode == "locked":
            frame_entry["rotor_angles"] = {"pm_rotor": rotor_angle_deg}
        timeline.append(frame_entry)

    circuits = [
        {
            "id": "stator_three_phase",
            "nodes": ["neutral", "phase_a", "phase_b", "phase_c"],
            "elements": [
                {"type": "resistor", "id": "R_a", "nodes": ["phase_a", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_a", "nodes": ["phase_a", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Va", "nodes": ["phase_a", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_a_pos", "inductor": "L_a", "region": "A_pos", "turns": turns_per_slot},
                {"type": "coil_link", "id": "coil_a_neg", "inductor": "L_a", "region": "A_neg", "turns": turns_per_slot},
                {"type": "resistor", "id": "R_b", "nodes": ["phase_b", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_b", "nodes": ["phase_b", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Vb", "nodes": ["phase_b", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_b_pos", "inductor": "L_b", "region": "B_pos", "turns": turns_per_slot},
                {"type": "coil_link", "id": "coil_b_neg", "inductor": "L_b", "region": "B_neg", "turns": turns_per_slot},
                {"type": "resistor", "id": "R_c", "nodes": ["phase_c", "neutral"], "resistance": phase_resistance},
                {"type": "inductor", "id": "L_c", "nodes": ["phase_c", "neutral"], "inductance": phase_inductance},
                {"type": "voltage_source", "id": "Vc", "nodes": ["phase_c", "neutral"], "value": 0.0},
                {"type": "coil_link", "id": "coil_c_pos", "inductor": "L_c", "region": "C_pos", "turns": turns_per_slot},
                {"type": "coil_link", "id": "coil_c_neg", "inductor": "L_c", "region": "C_neg", "turns": turns_per_slot},
            ],
        }
    ]

    initial_angle_deg = -load_angle_deg
    if mode == "spinup":
        initial_angle_deg = float(profile_cfg.get("spinup_initial_angle_deg", -load_angle_deg))

    scenario = {
        "version": "0.2",
        "units": "SI",
        "domain": {"Lx": domain_size, "Ly": domain_size, "nx": nx, "ny": ny},
        "boundary": {"type": "neumann"},
        "materials": [
            {"name": "air", "mu_r": 1.0},
            {"name": "stator_steel", "mu_r": 400.0},
            {"name": "rotor_steel", "mu_r": 800.0},
            {"name": "pm_magnet", "mu_r": magnet_mu_r},
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
                    "initial_angle_deg": initial_angle_deg,
                    "initial_speed_rpm": 0.0,
                    "torque_probe": torque_id,
                }
            ]
        },
        "timeline": timeline,
    }

    outputs = [
        {
            "type": "vtk_field_series",
            "id": f"{prefix}_series",
            "dir": "outputs",
            "basename": f"{prefix}_frame",
            "quantities": ["B", "H"],
        },
        {
            "type": "polyline_outlines",
            "id": f"{prefix}_outlines",
            "path": f"outputs/{prefix}_outlines.vtp",
        },
        {
            "type": "bore_avg_B",
            "id": f"{prefix}_bore",
            "path": f"outputs/{prefix}_bore.csv",
            "vertices": bore_probe_polygon,
        },
        {
            "type": "probe",
            "id": torque_id,
            "probe_type": "torque",
            "method": "stress_tensor",
            "loop": {"type": "polygon", "vertices": torque_probe_polygon},
            "path": f"outputs/{torque_id}.csv",
        },
    ]

    back_emf_outputs = []
    for phase, slots in phase_slot_map.items():
        pos_polygon = slots.get("pos")
        if not pos_polygon:
            continue
        back_emf_outputs.append(
            {
                "type": "back_emf_probe",
                "id": f"{prefix}_phase_{phase.lower()}_emf",
                "component": "Bmag",
                "region": {"type": "polygon", "vertices": pos_polygon},
                "path": f"outputs/{prefix}_phase_{phase.lower()}_emf.csv",
            }
        )

    outputs.extend(back_emf_outputs)

    if mode == "spinup":
        outputs.append(
            {
                "type": "mechanical_trace",
                "id": f"{prefix}_mechanical",
                "path": f"outputs/{prefix}_mechanical.csv",
                "rotors": ["pm_rotor"],
            }
        )

    scenario["outputs"] = outputs

    output_path.parent.mkdir(parents=True, exist_ok=True)
    quantized = quantize_numbers(scenario, ndigits=quantize_digits)
    output_path.write_text(json.dumps(quantized, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(DEFAULT_PROFILES.keys()), default="ci")
    parser.add_argument("--mode", choices=VALID_MODES, default="locked", help="Scenario mode: locked rotor or free spin-up")
    parser.add_argument("--cycles", type=int, help="Override the number of electrical cycles to emit")
    parser.add_argument("--frames-per-cycle", type=int, dest="frames_per_cycle", help="Override frames per electrical cycle")
    parser.add_argument("--out", type=Path, required=True, help="Output scenario JSON path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_scenario(
        args.profile,
        args.out,
        mode=args.mode,
        cycles_override=args.cycles,
        frames_per_cycle_override=args.frames_per_cycle,
    )


if __name__ == "__main__":
    main()
