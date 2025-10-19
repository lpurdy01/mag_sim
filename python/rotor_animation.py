#!/usr/bin/env python3
"""Reusable helpers for rendering rotor/stator animations."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import TwoSlopeNorm


@dataclass
class SlotGeometry:
    slot_id: str
    phase: str
    orientation: float
    turns: float
    vertices: np.ndarray


@dataclass
class RotorGeometry:
    pivot: np.ndarray
    polygons: List[np.ndarray]
    magnets: List[np.ndarray]


@dataclass
class CurrentSeries:
    times: np.ndarray
    effective: np.ndarray
    ampere_turns: np.ndarray


@dataclass
class ScenarioData:
    slots: List[SlotGeometry]
    rotor: RotorGeometry
    outer_polygons: List[np.ndarray]
    phase_series: Dict[str, CurrentSeries]


@dataclass
class MechanicalSeries:
    times: np.ndarray
    angles_rad: np.ndarray


def _as_array(points: Sequence[Sequence[float]]) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError("Expected Nx2 array for polygon vertices")
    return arr


def load_scenario(path: Path, rotor_name: str) -> ScenarioData:
    data = json.loads(path.read_text(encoding="utf-8"))
    regions = data.get("regions", [])
    rotors = data.get("rotors", [])
    sources = data.get("sources", [])
    magnets = data.get("magnet_regions", [])

    slot_geometries: List[SlotGeometry] = []
    for source in sources:
        if source.get("type") != "current_region":
            continue
        vertices = _as_array(source.get("vertices", []))
        slot_geometries.append(
            SlotGeometry(
                slot_id=str(source.get("id", "")),
                phase=str(source.get("phase", "")),
                orientation=float(source.get("orientation", 1.0)),
                turns=float(source.get("turns", 1.0)),
                vertices=vertices,
            )
        )

    rotor_entry = None
    for rotor in rotors:
        if rotor.get("name") == rotor_name:
            rotor_entry = rotor
            break
    if rotor_entry is None:
        raise ValueError(f"Rotor '{rotor_name}' not found in scenario {path}")

    pivot = np.asarray(rotor_entry.get("pivot", [0.0, 0.0]), dtype=float)

    rotor_polys: List[np.ndarray] = []
    for index in rotor_entry.get("polygons", []):
        if not isinstance(index, int):
            continue
        if index >= len(regions):
            continue
        region = regions[index]
        if "vertices" not in region:
            continue
        rotor_polys.append(_as_array(region["vertices"]))

    magnet_polys: List[np.ndarray] = []
    for index in rotor_entry.get("magnets", []):
        if not isinstance(index, int):
            continue
        if index >= len(magnets):
            continue
        magnet = magnets[index]
        if "vertices" in magnet:
            magnet_polys.append(_as_array(magnet["vertices"]))

    outer_polys: List[np.ndarray] = []
    for region in regions:
        if region.get("type") == "polygon" and "vertices" in region:
            outer_polys.append(_as_array(region["vertices"]))

    phase_series: Dict[str, CurrentSeries] = {}
    timeline = data.get("timeline", [])
    if isinstance(timeline, list) and timeline:
        samples: Dict[str, List[float]] = {}
        times: List[float] = []
        for entry in timeline:
            if not isinstance(entry, Mapping):
                continue
            time_value = float(entry.get("t", len(times)))
            if entry.get("phase_currents") and isinstance(entry["phase_currents"], Mapping):
                times.append(time_value)
                for phase, value in entry["phase_currents"].items():
                    samples.setdefault(str(phase), []).append(float(value))
        if times and samples:
            base_times = np.asarray(times, dtype=float)
            order = np.argsort(base_times)
            ordered_times = base_times[order]
            for phase, values in samples.items():
                if len(values) != len(ordered_times):
                    continue
                arr = np.asarray(values, dtype=float)[order]
                phase_series[phase] = CurrentSeries(
                    times=ordered_times,
                    effective=arr,
                    ampere_turns=arr,
                )

    scenario = ScenarioData(
        slots=slot_geometries,
        rotor=RotorGeometry(pivot=pivot, polygons=rotor_polys, magnets=magnet_polys),
        outer_polygons=outer_polys,
        phase_series=phase_series,
    )
    return scenario


def load_mechanical_trace(path: Path, rotor: str) -> MechanicalSeries:
    times: List[float] = []
    angles: List[float] = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("rotor") != rotor:
                continue
            try:
                time_val = float(row.get("time_s", row.get("time", 0.0)))
            except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
                raise ValueError(f"Invalid time entry in mechanical trace: {row!r}") from exc
            if "angle_rad" in row and row["angle_rad"] not in ("", None):
                angle = float(row["angle_rad"])
            else:
                angle_deg = float(row.get("angle_deg", 0.0))
                angle = math.radians(angle_deg)
            times.append(time_val)
            angles.append(angle)
    if not times:
        raise ValueError(f"Mechanical trace '{path}' did not contain samples for rotor '{rotor}'")
    arr_times = np.asarray(times, dtype=float)
    arr_angles = np.asarray(angles, dtype=float)
    order = np.argsort(arr_times)
    return MechanicalSeries(times=arr_times[order], angles_rad=arr_angles[order])


def load_circuit_trace(path: Path) -> Dict[str, CurrentSeries]:
    if not path.exists():
        raise FileNotFoundError(f"Circuit trace '{path}' not found")
    series: Dict[str, List[float]] = {}
    effective: Dict[str, List[float]] = {}
    ampere: Dict[str, List[float]] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Circuit trace '{path}' has no header")
        for row in reader:
            region = (row.get("region") or "").strip()
            coil = (row.get("coil") or "").strip()
            key = region or coil
            if not key:
                key = (row.get("circuit") or "").strip()
            if not key:
                continue
            time_val = float(row.get("time_s", row.get("time", 0.0)))
            eff = float(row.get("effective_current_A", row.get("current_A", 0.0)))
            amp = float(row.get("ampere_turns", eff))
            series.setdefault(key, []).append(time_val)
            effective.setdefault(key, []).append(eff)
            ampere.setdefault(key, []).append(amp)
    result: Dict[str, CurrentSeries] = {}
    for key, times in series.items():
        arr_times = np.asarray(times, dtype=float)
        order = np.argsort(arr_times)
        eff = np.asarray(effective[key], dtype=float)[order]
        amps = np.asarray(ampere[key], dtype=float)[order]
        result[key] = CurrentSeries(times=arr_times[order], effective=eff, ampere_turns=amps)
    return result


def extract_timeline_angles(scenario: ScenarioData, raw: Mapping[str, object], rotor: str) -> Optional[MechanicalSeries]:
    timeline = raw.get("timeline", [])
    if not isinstance(timeline, list) or not timeline:
        return None
    times: List[float] = []
    angles: List[float] = []
    for entry in timeline:
        if not isinstance(entry, Mapping):
            continue
        if "rotor_angles" not in entry:
            continue
        rotor_angles = entry["rotor_angles"]
        if not isinstance(rotor_angles, Mapping):
            continue
        if rotor not in rotor_angles:
            continue
        times.append(float(entry.get("t", len(times))))
        angles.append(math.radians(float(rotor_angles[rotor])))
    if not times:
        return None
    arr_times = np.asarray(times, dtype=float)
    order = np.argsort(arr_times)
    arr_angles = np.asarray(angles, dtype=float)[order]
    return MechanicalSeries(times=arr_times[order], angles_rad=arr_angles)


def rotate_polygon(points: np.ndarray, pivot: np.ndarray, angle_rad: float) -> np.ndarray:
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    translated = points - pivot
    rotated = np.empty_like(translated)
    rotated[:, 0] = translated[:, 0] * cos_a - translated[:, 1] * sin_a
    rotated[:, 1] = translated[:, 0] * sin_a + translated[:, 1] * cos_a
    return rotated + pivot


def _resample_series(series: CurrentSeries, sample_times: np.ndarray) -> np.ndarray:
    if series.times.size == 0:
        return np.zeros_like(sample_times)
    if series.times.size == 1:
        return np.full_like(sample_times, series.ampere_turns[0])
    return np.interp(sample_times, series.times, series.ampere_turns)


def build_slot_signal(
    slots: Sequence[SlotGeometry],
    sample_times: np.ndarray,
    circuit_series: Mapping[str, CurrentSeries],
    phase_series: Mapping[str, CurrentSeries],
) -> Dict[str, np.ndarray]:
    signals: Dict[str, np.ndarray] = {}
    for slot in slots:
        if slot.slot_id in circuit_series:
            series = circuit_series[slot.slot_id]
            values = _resample_series(series, sample_times)
        elif slot.phase and slot.phase in phase_series:
            series = phase_series[slot.phase]
            base = _resample_series(series, sample_times)
            orientation = slot.orientation if slot.orientation != 0.0 else 1.0
            values = base * slot.turns * orientation
        else:
            values = np.zeros_like(sample_times)
        signals[slot.slot_id] = values
    return signals


def compute_bounds(polygons: Iterable[np.ndarray]) -> Tuple[float, float, float, float]:
    min_x, max_x = math.inf, -math.inf
    min_y, max_y = math.inf, -math.inf
    for poly in polygons:
        if poly.size == 0:
            continue
        min_x = min(min_x, float(np.min(poly[:, 0])))
        max_x = max(max_x, float(np.max(poly[:, 0])))
        min_y = min(min_y, float(np.min(poly[:, 1])))
        max_y = max(max_y, float(np.max(poly[:, 1])))
    if not np.isfinite(min_x):
        min_x = max_x = min_y = max_y = 0.0
    padding_x = max(1e-3, 0.05 * (max_x - min_x + 1e-9))
    padding_y = max(1e-3, 0.05 * (max_y - min_y + 1e-9))
    return min_x - padding_x, max_x + padding_x, min_y - padding_y, max_y + padding_y


def render_animation(
    scenario: ScenarioData,
    mechanical: MechanicalSeries,
    slot_signals: Mapping[str, np.ndarray],
    *,
    sample_times: np.ndarray,
    gif_path: Optional[Path] = None,
    png_path: Optional[Path] = None,
    fps: int = 12,
    dpi: int = 200,
    max_frames: Optional[int] = None,
    figsize: Tuple[float, float] = (6.5, 6.5),
) -> None:
    n_frames = mechanical.times.size
    if max_frames is not None:
        n_frames = min(n_frames, max_frames)

    slot_values = list(slot_signals.values())
    max_ampere = max(float(np.max(np.abs(values))) for values in slot_values) if slot_values else 1.0
    if max_ampere <= 0.0:
        max_ampere = 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-max_ampere, vmax=max_ampere)
    cmap = plt.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")

    all_polys = [slot.vertices for slot in scenario.slots] + scenario.rotor.polygons + scenario.rotor.magnets + scenario.outer_polygons
    min_x, max_x, min_y, max_y = compute_bounds(all_polys)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    slot_patches: List[patches.Polygon] = []
    for slot in scenario.slots:
        patch = patches.Polygon(slot.vertices, closed=True, linewidth=0.4, edgecolor="#2f2f2f", facecolor=cmap(norm(0.0)))
        ax.add_patch(patch)
        slot_patches.append(patch)

    rotor_patches: List[patches.Polygon] = []
    for poly in scenario.rotor.polygons:
        patch = patches.Polygon(poly, closed=True, linewidth=0.6, edgecolor="#1b1b1b", facecolor="#a0a0a0", zorder=2)
        ax.add_patch(patch)
        rotor_patches.append(patch)

    magnet_patches: List[patches.Polygon] = []
    for poly in scenario.rotor.magnets:
        patch = patches.Polygon(poly, closed=True, linewidth=0.4, edgecolor="#333333", facecolor="#ffdd66", zorder=3)
        ax.add_patch(patch)
        magnet_patches.append(patch)

    time_text = ax.text(
        0.02,
        0.95,
        "",
        transform=ax.transAxes,
        fontsize=11,
        ha="left",
        va="top",
        color="#202020",
        bbox=dict(facecolor="white", alpha=0.7, boxstyle="round,pad=0.3"),
    )

    sample_indices = np.arange(n_frames)
    slot_arrays = [slot_signals[slot.slot_id][:n_frames] for slot in scenario.slots]

    def update(frame_index: int):
        angle = mechanical.angles_rad[min(frame_index, mechanical.angles_rad.size - 1)]
        time_val = mechanical.times[min(frame_index, mechanical.times.size - 1)]

        for patch, base in zip(rotor_patches, scenario.rotor.polygons):
            rotated = rotate_polygon(base, scenario.rotor.pivot, angle)
            patch.set_xy(rotated)
        for patch, base in zip(magnet_patches, scenario.rotor.magnets):
            rotated = rotate_polygon(base, scenario.rotor.pivot, angle)
            patch.set_xy(rotated)

        for patch, values in zip(slot_patches, slot_arrays):
            value = values[min(frame_index, values.size - 1)]
            patch.set_facecolor(cmap(norm(value)))

        rpm = angle * (60.0 / (2.0 * math.pi))
        time_text.set_text(f"t = {time_val:.4f} s\nangle = {math.degrees(angle):.2f}°\nω = {rpm:.1f} rpm")
        return slot_patches + rotor_patches + magnet_patches + [time_text]

    animation = FuncAnimation(fig, update, frames=sample_indices, interval=max(5, int(1000 / max(fps, 1))), blit=False)

    if gif_path:
        gif_path.parent.mkdir(parents=True, exist_ok=True)
        animation.save(gif_path, writer=PillowWriter(fps=fps))

    if png_path:
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, dpi=dpi, bbox_inches="tight")

    plt.close(fig)


__all__ = [
    "SlotGeometry",
    "RotorGeometry",
    "CurrentSeries",
    "ScenarioData",
    "MechanicalSeries",
    "load_scenario",
    "load_mechanical_trace",
    "load_circuit_trace",
    "extract_timeline_angles",
    "build_slot_signal",
    "render_animation",
]
