"""Generate three-phase induction motor scenarios with transient spin-up."""

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
        "cycles": 3,
        "slot_depth": 0.010,
        "slot_width_deg": 18.0,
        "bar_count": 6,
        "bar_width_deg": 12.0,
        "outer_segments": 24,
        "bore_segments": 18,
        "rotor_segments": 14,
        "torque_segments": 30,
        "bore_probe_segments": 24,
        "quantize_digits": 4,
    },
    "hires": {
        "nx": 401,
        "ny": 401,
        "frames_per_cycle": 180,
        "cycles": 6,
        "slot_depth": 0.010,
        "slot_width_deg": 12.0,
        "bar_count": 24,
        "bar_width_deg": 6.0,
        "outer_segments": 384,
        "bore_segments": 288,
        "rotor_segments": 256,
        "torque_segments": 384,
        "bore_probe_segments": 288,
        "quantize_digits": 6,
    },
}


VALID_MODES = ("spinup", "locked")


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


def quantize_numbers(value: object, ndigits: int = 6) -> object:
    """Round floats recursively so CI fixtures stay compact."""

    if isinstance(value, float):
        quantized = round(value, ndigits)
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


def compute_phase_currents(current_amp: float, omega: float, time: float) -> Tuple[float, float, float]:
    theta = omega * time
    ia = current_amp * math.sin(theta)
    ib = current_amp * math.sin(theta - 2.0 * math.pi / 3.0)
    ic = current_amp * math.sin(theta + 2.0 * math.pi / 3.0)
    return ia, ib, ic


def generate_scenario(
    profile: str,
    output_path: Path,
    *,
    mode: str = "spinup",
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

    domain_size = 0.14
    r_in = 0.040
    r_out = 0.055
    bore_fraction = 0.65

    current_amp = 42.0
    electrical_hz = 60.0
    omega = 2.0 * math.pi * electrical_hz

    rotor_inertia = 0.0011
    rotor_damping = 0.00008
    load_torque = 0.08

    if mode == "locked":
        cycles = max(1, cycles_override or 1)

    bar_count = int(profile_cfg.get("bar_count", 6))
    bar_width_deg = float(profile_cfg.get("bar_width_deg", 12.0))

    outer_segments = int(profile_cfg.get("outer_segments", 256))
    bore_segments = int(profile_cfg.get("bore_segments", 192))
    rotor_segments = int(profile_cfg.get("rotor_segments", 160))
    torque_segments = int(profile_cfg.get("torque_segments", max(64, rotor_segments)))
    bore_probe_segments = int(profile_cfg.get("bore_probe_segments", max(64, bore_segments)))
    quantize_digits = int(profile_cfg.get("quantize_digits", 6))

    total_frames = frames_per_cycle * cycles
    dt = 1.0 / electrical_hz / frames_per_cycle

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
    for _phase, _orientation, angle, _slot_id in slot_definitions:
        slot_polygons.append(build_slot_polygon(angle, slot_inner_radius, slot_outer_radius, slot_width_deg))

    regions: List[Dict[str, object]] = [
        {"type": "uniform", "material": "air"},
        {"type": "polygon", "material": "stator_steel", "vertices": outer_polygon},
    ]
    next_polygon_index = 1
    for polygon in slot_polygons:
        regions.append({"type": "polygon", "material": "air", "vertices": polygon})
        next_polygon_index += 1
    regions.append({"type": "polygon", "material": "air", "vertices": bore_polygon})
    next_polygon_index += 1

    rotor_radius = r_in * 0.58
    rotor_polygon = build_circle(rotor_radius, rotor_segments)
    rotor_polygon_index = next_polygon_index
    regions.append({"type": "polygon", "material": "rotor_core", "vertices": rotor_polygon})
    next_polygon_index += 1

    bar_inner_radius = rotor_radius * 0.35
    bar_outer_radius = rotor_radius * 0.98
    bar_polygons: List[List[List[float]]] = []
    for bar_idx in range(bar_count):
        angle = 2.0 * math.pi * bar_idx / bar_count
        bar_polygons.append(build_slot_polygon(angle, bar_inner_radius, bar_outer_radius, bar_width_deg))

    rotor_polygon_indices = [rotor_polygon_index]
    for polygon in bar_polygons:
        regions.append({"type": "polygon", "material": "rotor_bar", "vertices": polygon})
        rotor_polygon_indices.append(next_polygon_index)
        next_polygon_index += 1

    magnet_regions: List[Dict[str, object]] = []

    rotors = [
        {
            "name": "induction_rotor",
            "pivot": [0.0, 0.0],
            "polygons": rotor_polygon_indices,
        }
    ]

    phase_turns = 100.0
    turns_per_slot = phase_turns / 2.0
    slot_fill_fraction = 0.55

    sources = []
    phase_slot_map: Dict[str, Dict[str, List[List[float]]]] = {"A": {}, "B": {}, "C": {}}
    for phase, orientation, angle, slot_id in slot_definitions:
        polygon = build_slot_polygon(angle, slot_inner_radius, slot_outer_radius, slot_width_deg)
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

    materials = [
        {"name": "air", "mu_r": 1.0, "sigma": 0.0},
        {"name": "stator_steel", "mu_r": 600.0, "sigma": 0.0},
        {"name": "rotor_core", "mu_r": 400.0, "sigma": 5.0e5},
        {"name": "rotor_bar", "mu_r": 1.0, "sigma": 1.8e7},
    ]

    timeline = []
    for frame_idx in range(total_frames):
        time = frame_idx * dt
        ia, ib, ic = compute_phase_currents(current_amp, omega, time)
        entry = {
            "t": time,
            "phase_currents": {
                "A": ia,
                "B": ib,
                "C": ic,
            },
        }
        if mode == "locked":
            electrical_angle_deg = math.degrees(omega * time)
            entry["rotor_angles"] = {"induction_rotor": electrical_angle_deg}
        timeline.append(entry)

    mechanical = {
        "rotors": [
            {
                "name": "induction_rotor",
                "inertia": rotor_inertia,
                "damping": rotor_damping,
                "load_torque": load_torque,
                "initial_angle_deg": -20.0 if mode == "spinup" else 0.0,
                "initial_speed_rpm": 0.0,
                "torque_probe": "induction_motor_torque",
            }
        ]
    }

    scenario = {
        "version": "0.2",
        "units": "SI",
        "domain": {"Lx": domain_size, "Ly": domain_size, "nx": nx, "ny": ny},
        "boundary": {"type": "neumann"},
        "materials": materials,
        "regions": regions,
        "magnet_regions": magnet_regions,
        "rotors": rotors,
        "sources": sources,
        "mechanical": mechanical,
        "transient": {"dt": dt, "n_steps": total_frames},
        "timeline": timeline,
    }

    outputs = [
        {
            "type": "vtk_field_series",
            "id": "induction_motor_series",
            "dir": "outputs",
            "basename": "induction_motor_frame",
            "quantities": ["B", "H"],
        },
        {
            "type": "polyline_outlines",
            "id": "induction_motor_outlines",
            "path": "outputs/induction_motor_outlines.vtp",
        },
        {
            "type": "bore_avg_B",
            "id": "induction_motor_bore",
            "path": "outputs/induction_motor_bore.csv",
            "vertices": build_circle(r_in * bore_fraction, bore_probe_segments),
        },
        {
            "type": "probe",
            "id": "induction_motor_torque",
            "probe_type": "torque",
            "method": "stress_tensor",
            "loop": {"type": "polygon", "vertices": build_circle(rotor_radius + 0.003, torque_segments)},
            "path": "outputs/induction_motor_torque.csv",
        },
    ]

    for phase, slots in phase_slot_map.items():
        pos_polygon = slots.get("pos")
        if not pos_polygon:
            continue
        outputs.append(
            {
                "type": "back_emf_probe",
                "id": f"induction_motor_phase_{phase.lower()}_emf",
                "component": "Bmag",
                "region": {"type": "polygon", "vertices": pos_polygon},
                "path": f"outputs/induction_motor_phase_{phase.lower()}_emf.csv",
            }
        )

    outputs.append(
        {
            "type": "mechanical_trace",
            "id": "induction_motor_mechanical",
            "path": "outputs/induction_motor_mechanical.csv",
            "rotors": ["induction_rotor"],
        }
    )

    scenario["outputs"] = outputs

    output_path.parent.mkdir(parents=True, exist_ok=True)
    quantized = quantize_numbers(scenario, ndigits=quantize_digits)
    output_path.write_text(json.dumps(quantized, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(DEFAULT_PROFILES.keys()), default="ci")
    parser.add_argument("--mode", choices=VALID_MODES, default="spinup", help="Scenario mode: free spin-up or locked rotor")
    parser.add_argument("--cycles", type=int, help="Override the number of electrical cycles to emit")
    parser.add_argument(
        "--frames-per-cycle", type=int, dest="frames_per_cycle", help="Override frames per electrical cycle"
    )
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
