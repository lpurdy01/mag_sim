"""Create animations for motor field maps with current and rotor overlays."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import vtk  # type: ignore
from matplotlib import patheffects
from matplotlib.colors import LogNorm
from matplotlib.patches import FancyArrowPatch, Polygon as MplPolygon
from vtk.util.numpy_support import vtk_to_numpy  # type: ignore

from rotor_animation import (  # type: ignore
    CurrentSeries as RotorCurrentSeries,
    MechanicalSeries,
    load_circuit_trace,
    load_mechanical_trace,
)


@dataclass
class CurrentSeriesData:
    label: str
    times: np.ndarray
    values: np.ndarray


@dataclass
class RotorPolygon:
    vertices: np.ndarray
    material: str


@dataclass
class RotorOverlay:
    name: str
    pivot: np.ndarray
    rotor_polygons: List[RotorPolygon]
    magnet_polygons: List[RotorPolygon]


def resolve_optional_path(base: Path, candidate: str) -> Path:
    path = Path(candidate)
    if path.is_absolute():
        return path
    search_roots: Sequence[Path] = (Path.cwd(), base, base.parent)
    for root in search_roots:
        trial = (root / path).resolve()
        if trial.exists():
            return trial
    return (Path.cwd() / path).resolve()


def rotate_polygon(points: np.ndarray, pivot: np.ndarray, angle_rad: float) -> np.ndarray:
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    translated = points - pivot
    rotated = np.empty_like(translated)
    rotated[:, 0] = translated[:, 0] * cos_a - translated[:, 1] * sin_a
    rotated[:, 1] = translated[:, 0] * sin_a + translated[:, 1] * cos_a
    return rotated + pivot


def rotor_patch_style(material: str, *, is_magnet: bool = False) -> Tuple[Tuple[float, float, float, float], str, float]:
    key = material.lower()
    if is_magnet:
        return (1.0, 0.74, 0.3, 0.72), "#7a4200", 1.1
    if "bar" in key or "conductor" in key:
        return (0.96, 0.73, 0.32, 0.75), "#915b0f", 1.05
    if "air" in key:
        return (0.75, 0.75, 0.75, 0.35), "#404040", 0.8
    return (0.86, 0.86, 0.86, 0.68), "#202020", 1.0


def extract_rotor_overlay(scenario: Dict[str, object], rotor_name: Optional[str]) -> Optional[RotorOverlay]:
    rotors = scenario.get("rotors", [])
    if not isinstance(rotors, list) or not rotors:
        return None
    selected: Optional[Dict[str, object]] = None
    if rotor_name is not None:
        for entry in rotors:
            if isinstance(entry, dict) and entry.get("name") == rotor_name:
                selected = entry
                break
    if selected is None:
        for entry in rotors:
            if isinstance(entry, dict):
                selected = entry
                break
    if selected is None:
        return None

    name = str(selected.get("name", rotor_name or "rotor"))
    pivot = np.asarray(selected.get("pivot", [0.0, 0.0]), dtype=float)

    regions = scenario.get("regions", [])
    rotor_polys: List[RotorPolygon] = []
    for index in selected.get("polygons", []) or []:
        if not isinstance(index, int):
            continue
        if not isinstance(regions, list) or index >= len(regions):
            continue
        region = regions[index]
        if not isinstance(region, dict):
            continue
        vertices = region.get("vertices")
        if isinstance(vertices, list) and vertices:
            rotor_polys.append(
                RotorPolygon(vertices=np.asarray(vertices, dtype=float), material=str(region.get("material", "")))
            )

    magnet_defs = scenario.get("magnet_regions", [])
    magnet_polys: List[RotorPolygon] = []
    for index in selected.get("magnets", []) or []:
        if not isinstance(index, int):
            continue
        if not isinstance(magnet_defs, list) or index >= len(magnet_defs):
            continue
        magnet = magnet_defs[index]
        if not isinstance(magnet, dict):
            continue
        vertices = magnet.get("vertices")
        if isinstance(vertices, list) and vertices:
            magnet_polys.append(
                RotorPolygon(vertices=np.asarray(vertices, dtype=float), material=str(magnet.get("material", "")))
            )

    return RotorOverlay(name=name, pivot=pivot, rotor_polygons=rotor_polys, magnet_polygons=magnet_polys)


def extract_timeline_mechanical(scenario: Dict[str, object], rotor_name: str) -> Optional[MechanicalSeries]:
    timeline = scenario.get("timeline", [])
    if not isinstance(timeline, list) or not timeline:
        return None
    times: List[float] = []
    angles: List[float] = []
    for entry in timeline:
        if not isinstance(entry, dict):
            continue
        rotor_angles = entry.get("rotor_angles")
        if not isinstance(rotor_angles, dict):
            continue
        if rotor_name not in rotor_angles:
            continue
        times.append(float(entry.get("t", entry.get("time", len(times)))))
        angles.append(math.radians(float(rotor_angles[rotor_name])))
    if not times:
        return None
    arr_times = np.asarray(times, dtype=float)
    order = np.argsort(arr_times)
    arr_angles = np.asarray(angles, dtype=float)[order]
    return MechanicalSeries(times=arr_times[order], angles_rad=arr_angles)


def infer_output_path(scenario: Dict[str, object], output_type: str, scenario_path: Path) -> Optional[Path]:
    outputs = scenario.get("outputs", [])
    if not isinstance(outputs, list):
        return None
    for entry in outputs:
        if not isinstance(entry, dict):
            continue
        if entry.get("type") != output_type:
            continue
        raw_path = entry.get("path") or entry.get("file")
        if isinstance(raw_path, str) and raw_path:
            return resolve_optional_path(scenario_path.parent, raw_path)
    return None


def load_current_series(
    scenario: Dict[str, object],
    circuit_trace: Optional[Path],
) -> List[CurrentSeriesData]:
    timeline = scenario.get("timeline", [])
    series_list: List[CurrentSeriesData] = []
    if isinstance(timeline, list) and timeline:
        records: List[Tuple[float, Dict[str, float]]] = []
        labels: List[str] = []
        for entry in timeline:
            if not isinstance(entry, dict):
                continue
            phases = entry.get("phase_currents")
            if not isinstance(phases, dict):
                continue
            time = float(entry.get("t", entry.get("time", len(records))))
            casted = {str(key): float(value) for key, value in phases.items()}
            records.append((time, casted))
            for key in casted.keys():
                if key not in labels:
                    labels.append(key)
        if records and labels:
            times = np.asarray([item[0] for item in records], dtype=float)
            order = np.argsort(times)
            ordered_times = times[order]
            for label in labels:
                values = np.zeros_like(ordered_times)
                for idx, (_, phases) in enumerate(records):
                    if label in phases:
                        values[idx] = phases[label]
                ordered_values = values[order]
                series_list.append(CurrentSeriesData(label=label, times=ordered_times, values=ordered_values))
    if series_list:
        return series_list

    if circuit_trace is None or not circuit_trace.exists():
        return []
    try:
        trace = load_circuit_trace(circuit_trace)
    except Exception as exc:  # pragma: no cover - surfaced to CLI
        raise RuntimeError(f"Failed to load circuit trace '{circuit_trace}': {exc}") from exc
    ordered = sorted(trace.items(), key=lambda item: item[0])
    for label, series in ordered:
        if not isinstance(series, RotorCurrentSeries):
            continue
        series_list.append(
            CurrentSeriesData(
                label=str(label) or "current",
                times=np.asarray(series.times, dtype=float),
                values=np.asarray(series.ampere_turns, dtype=float),
            )
        )
    return series_list


def sample_rotor_angles(
    mechanical: Optional[MechanicalSeries], sample_times: np.ndarray
) -> Optional[np.ndarray]:
    if mechanical is None or mechanical.times.size == 0:
        return None
    base_times = mechanical.times
    base_angles = mechanical.angles_rad
    if base_times.size == 1:
        return np.full_like(sample_times, base_angles[0], dtype=float)
    return np.interp(
        sample_times,
        base_times,
        base_angles,
        left=base_angles[0],
        right=base_angles[-1],
    )


def point_in_polygon(x: float, y: float, vertices: np.ndarray) -> bool:
    inside = False
    count = vertices.shape[0]
    for i in range(count):
        j = (i - 1) % count
        xi, yi = vertices[i]
        xj, yj = vertices[j]
        intersect = ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        if intersect:
            inside = not inside
    return inside


def polygon_area(vertices: np.ndarray) -> float:
    if vertices.shape[0] < 3:
        return 0.0
    xs = vertices[:, 0]
    ys = vertices[:, 1]
    return 0.5 * float(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1)))


def polygon_centroid(vertices: np.ndarray) -> Tuple[float, float]:
    if vertices.shape[0] < 3:
        return float(np.mean(vertices[:, 0])), float(np.mean(vertices[:, 1]))
    area = polygon_area(vertices)
    if abs(area) < 1e-12:
        return float(np.mean(vertices[:, 0])), float(np.mean(vertices[:, 1]))
    xs = vertices[:, 0]
    ys = vertices[:, 1]
    cross = xs * np.roll(ys, -1) - ys * np.roll(xs, -1)
    factor = 1.0 / (6.0 * area)
    cx = float(np.sum((xs + np.roll(xs, -1)) * cross) * factor)
    cy = float(np.sum((ys + np.roll(ys, -1)) * cross) * factor)
    return cx, cy


def extract_outline_polygons(scenario: Dict[str, object]) -> List[Dict[str, np.ndarray]]:
    regions = scenario.get("regions", [])
    polygons: List[Dict[str, object]] = []
    for idx, region in enumerate(regions):
        if not isinstance(region, dict):
            continue
        if region.get("type") != "polygon":
            continue
        vertices = region.get("vertices")
        if not isinstance(vertices, list) or len(vertices) < 3:
            continue
        arr = np.asarray(vertices, dtype=float)
        polygons.append(
            {
                "index": idx,
                "material": region.get("material", ""),
                "vertices": arr,
                "area": abs(polygon_area(arr)),
            }
        )

    if not polygons:
        return []

    outlines: List[Dict[str, np.ndarray]] = []
    used_indices: set[int] = set()

    def add_outline(category: str, poly: Dict[str, object]) -> None:
        if poly["index"] in used_indices:
            return
        outlines.append({"category": category, "vertices": np.asarray(poly["vertices"], dtype=float)})
        used_indices.add(poly["index"])

    stator_candidates = [poly for poly in polygons if poly.get("material") not in {"", "air"}]
    if stator_candidates:
        stator_poly = max(stator_candidates, key=lambda poly: float(poly["area"]))
    else:
        stator_poly = max(polygons, key=lambda poly: float(poly["area"]))
    add_outline("stator", stator_poly)

    air_polys = [poly for poly in polygons if poly.get("material") == "air"]
    bore_poly = None
    if air_polys:
        bore_poly = max(air_polys, key=lambda poly: float(poly["area"]))
        add_outline("bore", bore_poly)

    for poly in polygons:
        if poly["index"] in used_indices:
            continue
        add_outline("slot", poly)

    return outlines


def extract_slot_annotations(scenario: Dict[str, object]) -> List[Dict[str, object]]:
    slots: List[Dict[str, object]] = []
    for source in scenario.get("sources", []):
        if not isinstance(source, dict):
            continue
        if source.get("type") not in {"current_region", "current-region"}:
            continue
        vertices = source.get("vertices")
        if not isinstance(vertices, list) or len(vertices) < 3:
            continue
        arr = np.asarray(vertices, dtype=float)
        phase = str(source.get("phase", "")) if source.get("phase") is not None else ""
        label = source.get("label") or source.get("id") or phase
        orientation = float(source.get("orientation", 1.0))
        if not source.get("label") and phase:
            sign = "+" if orientation >= 0.0 else "-"
            label = f"{phase}{sign}"
        slots.append({"phase": phase, "label": label, "vertices": arr})
    return slots


def load_scenario(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_pvd(path: Path) -> List[Tuple[float, Path]]:
    root = ET.fromstring(path.read_text(encoding="utf-8"))
    datasets: List[Tuple[float, Path]] = []
    base_dir = path.parent
    for dataset in root.findall(".//DataSet"):
        timestep = float(dataset.attrib.get("timestep", "0.0"))
        file_attr = dataset.attrib.get("file")
        if not file_attr:
            continue
        datasets.append((timestep, (base_dir / file_attr).resolve()))
    datasets.sort(key=lambda item: item[0])
    return datasets


def extract_bore_polygon(scenario: Dict[str, object]) -> np.ndarray:
    outputs = scenario.get("outputs", [])
    for output in outputs:
        if isinstance(output, dict) and output.get("type") in {"bore_avg_B", "bore_avg_b"}:
            vertices = output.get("vertices") or output.get("polygon")
            if isinstance(vertices, dict):
                vertices = vertices.get("vertices")
            if isinstance(vertices, list) and vertices:
                return np.asarray(vertices, dtype=float)
    # Fallback: use bore radius from geometry
    domain = scenario.get("domain", {})
    nx = int(domain.get("nx", 1))
    ny = int(domain.get("ny", 1))
    Lx = float(domain.get("Lx", 1.0))
    dx = Lx / max(nx - 1, 1)
    radius = 0.5 * Lx * 0.3
    segments = 64
    vertices = np.array(
        [[radius * math.cos(2.0 * math.pi * i / segments), radius * math.sin(2.0 * math.pi * i / segments)] for i in range(segments)],
        dtype=float,
    )
    return vertices


def compute_bore_series(
    datasets: List[Tuple[float, Path]], bore_polygon: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    times = []
    bx_values = []
    by_values = []
    for time, vti_path in datasets:
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(str(vti_path))
        reader.Update()
        data = reader.GetOutput()
        cell_data = data.GetCellData()
        bx_array = cell_data.GetArray("Bx")
        by_array = cell_data.GetArray("By")
        if bx_array is None or by_array is None:
            raise RuntimeError(f"VTK file {vti_path} is missing Bx/By arrays")
        bx = vtk_to_numpy(bx_array)
        by = vtk_to_numpy(by_array)
        dims = data.GetDimensions()
        spacing = data.GetSpacing()
        origin = data.GetOrigin()
        cell_nx = max(dims[0] - 1, 1)
        cell_ny = max(dims[1] - 1, 1)
        bx = bx.reshape((cell_ny, cell_nx))
        by = by.reshape((cell_ny, cell_nx))
        mask = np.zeros((cell_ny, cell_nx), dtype=bool)
        for j in range(cell_ny):
            y = origin[1] + (j + 0.5) * spacing[1]
            for i in range(cell_nx):
                x = origin[0] + (i + 0.5) * spacing[0]
                if point_in_polygon(x, y, bore_polygon):
                    mask[j, i] = True
        if not np.any(mask):
            raise RuntimeError(f"Bore polygon did not overlap any cells in {vti_path}")
        avg_bx = float(np.mean(bx[mask]))
        avg_by = float(np.mean(by[mask]))
        times.append(time)
        bx_values.append(avg_bx)
        by_values.append(avg_by)
    return np.asarray(times, dtype=float), np.asarray(bx_values, dtype=float), np.asarray(by_values, dtype=float)


def load_field_frames(
    datasets: List[Tuple[float, Path]]
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, np.ndarray]]]:
    frames: List[Dict[str, np.ndarray]] = []
    x_coords: Optional[np.ndarray] = None
    y_coords: Optional[np.ndarray] = None
    for _, vti_path in datasets:
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(str(vti_path))
        reader.Update()
        data = reader.GetOutput()
        cell_data = data.GetCellData()
        bx_array = cell_data.GetArray("Bx")
        by_array = cell_data.GetArray("By")
        if bx_array is None or by_array is None:
            raise RuntimeError(f"VTK file {vti_path} is missing Bx/By arrays")
        bx = vtk_to_numpy(bx_array)
        by = vtk_to_numpy(by_array)
        dims = data.GetDimensions()
        spacing = data.GetSpacing()
        origin = data.GetOrigin()
        cell_nx = max(dims[0] - 1, 1)
        cell_ny = max(dims[1] - 1, 1)
        bx = bx.reshape((cell_ny, cell_nx))
        by = by.reshape((cell_ny, cell_nx))
        mag = np.hypot(bx, by)
        if x_coords is None or y_coords is None:
            x_centers = origin[0] + (np.arange(cell_nx) + 0.5) * spacing[0]
            y_centers = origin[1] + (np.arange(cell_ny) + 0.5) * spacing[1]
            x_coords, y_coords = np.meshgrid(x_centers, y_centers)
        frames.append({"bx": bx, "by": by, "mag": mag})
    if x_coords is None or y_coords is None:
        raise RuntimeError("No field frames could be loaded")
    return x_coords, y_coords, frames


def build_animation(
    save_path: Path,
    html_path: Optional[Path],
    still_path: Optional[Path],
    scenario: Dict[str, object],
    datasets: List[Tuple[float, Path]],
    current_series: List[CurrentSeriesData],
    bore_times: np.ndarray,
    bore_bx: np.ndarray,
    bore_by: np.ndarray,
    bore_polygon: np.ndarray,
    rotor_overlay: Optional[RotorOverlay],
    rotor_angles: Optional[np.ndarray],
    fps: int,
    width: int,
    log_scale: bool,
) -> None:
    x_coords, y_coords, field_frames = load_field_frames(datasets)
    bore_magnitude = np.hypot(bore_bx, bore_by)
    max_mag = float(np.max(bore_magnitude)) if bore_magnitude.size > 0 else 1.0
    bore_angles = np.arctan2(bore_by, bore_bx)
    outlines = extract_outline_polygons(scenario)
    slot_annotations = extract_slot_annotations(scenario)
    pivot_point = rotor_overlay.pivot if rotor_overlay is not None else np.array([0.0, 0.0])
    stator_radius = None
    for outline in outlines:
        if outline.get("category") != "stator":
            continue
        vertices = outline.get("vertices")
        if isinstance(vertices, np.ndarray) and vertices.size:
            offsets = vertices - pivot_point
            stator_radius = float(np.max(np.linalg.norm(offsets, axis=1)))
            break

    fig_width_in = width / 100.0
    fig_height_in = fig_width_in * 0.9
    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    gs = fig.add_gridspec(3, 2, height_ratios=[6, 2, 2], width_ratios=[20, 1], hspace=0.3, wspace=0.05)

    ax_field = fig.add_subplot(gs[0, 0])
    ax_field.set_aspect("equal")
    extent = [x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()]
    all_mags = [frame["mag"] for frame in field_frames]
    mag_max = max(float(np.max(mag)) for mag in all_mags)
    mag_min = min(float(np.min(mag)) for mag in all_mags)
    if log_scale:
        positive_samples = [mag[mag > 0.0] for mag in all_mags]
        positives = (
            np.concatenate([sample for sample in positive_samples if sample.size])
            if any(sample.size for sample in positive_samples)
            else np.array([], dtype=float)
        )
        if positives.size:
            vmin = float(np.percentile(positives, 5))
            vmax = float(np.max(positives))
        else:
            vmin, vmax = 1e-6, 1.0
        norm = LogNorm(vmin=max(vmin, 1e-9), vmax=max(vmax, vmin * 10.0))
        field_to_display = lambda mag: np.maximum(mag, norm.vmin)  # noqa: E731
    else:
        norm = None
        field_to_display = lambda mag: mag  # noqa: E731
    first_mag = all_mags[0]
    im = ax_field.imshow(
        field_to_display(first_mag),
        origin="lower",
        extent=extent,
        cmap="viridis",
        norm=norm,
        vmin=None if log_scale else mag_min,
        vmax=None if log_scale else mag_max,
    )
    ax_field.set_xlabel("x [m]")
    ax_field.set_ylabel("y [m]")

    quiver_step = max(1, int(max(first_mag.shape) / 32))
    sampled_x = x_coords[::quiver_step, ::quiver_step]
    sampled_y = y_coords[::quiver_step, ::quiver_step]
    all_vectors = [np.hypot(frame["bx"], frame["by"]) for frame in field_frames]
    vector_max = max(float(np.max(mag)) for mag in all_vectors)
    dx = float(np.mean(np.diff(sampled_x[0, :]))) if sampled_x.shape[1] > 1 else (extent[1] - extent[0])
    dy = float(np.mean(np.diff(sampled_y[:, 0]))) if sampled_y.shape[0] > 1 else (extent[3] - extent[2])
    base_spacing = max(dx, dy, 1e-6)
    target_arrow = base_spacing * 0.75
    quiver_scale = vector_max / target_arrow if vector_max > 0.0 else 1.0
    quiver = ax_field.quiver(
        sampled_x,
        sampled_y,
        field_frames[0]["bx"][::quiver_step, ::quiver_step],
        field_frames[0]["by"][::quiver_step, ::quiver_step],
        color="white",
        scale=quiver_scale,
        scale_units="xy",
        angles="xy",
        width=0.004,
        pivot="mid",
    )

    rotor_patches: List[MplPolygon] = []
    magnet_patches: List[MplPolygon] = []
    rotor_base_polys: List[np.ndarray] = []
    magnet_base_polys: List[np.ndarray] = []
    if rotor_overlay is not None:
        for poly in rotor_overlay.rotor_polygons:
            base = np.asarray(poly.vertices, dtype=float)
            rotor_base_polys.append(base)
            face, edge, lw = rotor_patch_style(poly.material)
            patch = MplPolygon(
                base,
                closed=True,
                facecolor=face,
                edgecolor=edge,
                linewidth=lw,
                zorder=6,
            )
            ax_field.add_patch(patch)
            rotor_patches.append(patch)
        for poly in rotor_overlay.magnet_polygons:
            base = np.asarray(poly.vertices, dtype=float)
            magnet_base_polys.append(base)
            face, edge, lw = rotor_patch_style(poly.material, is_magnet=True)
            patch = MplPolygon(
                base,
                closed=True,
                facecolor=face,
                edgecolor=edge,
                linewidth=lw,
                zorder=7,
            )
            ax_field.add_patch(patch)
            magnet_patches.append(patch)

    outline_styles = {
        "stator": {"color": "white", "linewidth": 2.0},
        "bore": {"color": "tab:orange", "linewidth": 1.8},
        "slot": {"color": "white", "linewidth": 1.2},
    }
    for outline in outlines:
        vertices = outline.get("vertices")
        if not isinstance(vertices, np.ndarray) or vertices.size == 0:
            continue
        closed = np.vstack([vertices, vertices[0]]) if vertices.shape[0] >= 2 else vertices
        style = outline_styles.get(outline.get("category", "slot"), {"color": "white", "linewidth": 1.0})
        line, = ax_field.plot(
            closed[:, 0],
            closed[:, 1],
            color=style["color"],
            lw=style["linewidth"],
            zorder=5,
        )
        line.set_path_effects(
            [
                patheffects.Stroke(linewidth=style["linewidth"] + 1.0, foreground="black"),
                patheffects.Normal(),
            ]
        )

    phase_colors = {
        "A": "tab:red",
        "B": "tab:green",
        "C": "tab:blue",
        "field": "tab:purple",
        "armature": "tab:orange",
    }
    label_ring = stator_radius if stator_radius is not None else (np.max(np.linalg.norm(bore_polygon, axis=1)) if bore_polygon.size else 0.1)
    label_ring = float(label_ring)
    label_offset = max(label_ring * 0.12, 0.01)
    for slot in slot_annotations:
        vertices = slot.get("vertices")
        if not isinstance(vertices, np.ndarray) or vertices.size == 0:
            continue
        cx, cy = polygon_centroid(vertices)
        label = str(slot.get("label", ""))
        color = phase_colors.get(str(slot.get("phase", "")), "white")
        direction = np.array([cx, cy]) - pivot_point
        dist = float(np.linalg.norm(direction))
        if dist < 1e-9:
            direction = np.array([1.0, 0.0])
            dist = 1e-9
        else:
            direction = direction / dist
        text_radius = max(dist + label_offset, label_ring + label_offset)
        target = pivot_point + direction * text_radius
        ax_field.annotate(
            label,
            xy=(cx, cy),
            xytext=(target[0], target[1]),
            textcoords="data",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color=color,
            arrowprops={"arrowstyle": "-", "color": color, "lw": 1.1},
            bbox={
                "boxstyle": "round,pad=0.28",
                "facecolor": (0.05, 0.05, 0.05, 0.78),
                "edgecolor": color,
                "linewidth": 0.8,
            },
            zorder=8,
        )

    ax_field.set_title("Magnetic flux density magnitude")
    cbar = fig.colorbar(im, cax=fig.add_subplot(gs[0, 1]))
    cbar.set_label("|B| [T]")

    fig.add_subplot(gs[1, 1]).axis("off")

    ax_compass = fig.add_subplot(gs[1, 0])
    bore_radius = np.max(np.linalg.norm(bore_polygon, axis=1)) if bore_polygon.size else 1.0
    circle = plt.Circle((0.0, 0.0), bore_radius, fill=False, color="tab:gray", lw=1.2)
    ax_compass.add_patch(circle)
    ax_compass.set_aspect("equal", adjustable="datalim")
    padding = bore_radius * 1.2
    ax_compass.set_xlim(-padding, padding)
    ax_compass.set_ylim(-padding, padding)
    ax_compass.axis("off")
    arrow = FancyArrowPatch((0.0, 0.0), (0.0, 0.0), color="tab:orange", arrowstyle="->", linewidth=2.0)
    ax_compass.add_patch(arrow)
    angle_text = ax_compass.text(0.02, 0.9, "", transform=ax_compass.transAxes, ha="left", va="top")

    ax_currents = fig.add_subplot(gs[2, 0])
    fig.add_subplot(gs[2, 1]).axis("off")

    lines: List[object] = []
    phase_colors = {"A": "tab:red", "B": "tab:green", "C": "tab:blue"}
    for idx, series in enumerate(current_series):
        if series.times.size == 0:
            continue
        color = phase_colors.get(series.label.upper(), f"C{idx}")
        line, = ax_currents.plot(series.times, series.values, label=series.label, color=color)
        lines.append(line)
    ax_currents.set_xlabel("Time [s]")
    ax_currents.set_ylabel("Current [A-turns]")
    if lines:
        ax_currents.legend(loc="upper right")
    ax_currents.grid(True, alpha=0.3)
    if current_series and any(series.times.size for series in current_series):
        min_time = min(float(np.min(series.times)) for series in current_series if series.times.size)
        max_time = max(float(np.max(series.times)) for series in current_series if series.times.size)
    elif bore_times.size:
        min_time = float(bore_times[0])
        max_time = float(bore_times[-1])
    else:
        min_time = 0.0
        max_time = 1.0
    if bore_times.size:
        frame_start = float(bore_times[0])
        frame_end = float(bore_times[-1])
        min_time = min(min_time, frame_start)
        max_time = min(max_time, frame_end)
    if math.isclose(min_time, max_time):
        max_time = min_time + 1e-6
    ax_currents.set_xlim(min_time, max_time)
    marker_initial = bore_times[0] if bore_times.size else min_time
    marker_line = ax_currents.axvline(marker_initial, color="k", linestyle="--")

    def update(frame_idx: int) -> List[object]:
        field = field_frames[frame_idx]
        im.set_data(field_to_display(field["mag"]))
        quiver.set_UVC(
            field["bx"][::quiver_step, ::quiver_step],
            field["by"][::quiver_step, ::quiver_step],
        )
        time = bore_times[frame_idx]
        magnitude = bore_magnitude[frame_idx]
        angle = bore_angles[frame_idx]
        scale = 0.8 * bore_radius * (magnitude / max_mag if max_mag > 0 else 1.0)
        dx = scale * math.cos(angle)
        dy = scale * math.sin(angle)
        arrow.set_positions((0.0, 0.0), (dx, dy))
        angle_deg = math.degrees(angle)
        angle_text.set_text(f"t={time:.4f} s\n|B|={magnitude:.3f} T\nangle={angle_deg:.1f}°")
        marker_line.set_xdata([time, time])
        artists: List[object] = [im, quiver, arrow, angle_text, marker_line]
        if rotor_overlay is not None and rotor_angles is not None and rotor_patches:
            rotor_angle = rotor_angles[min(frame_idx, rotor_angles.size - 1)]
            for patch, base in zip(rotor_patches, rotor_base_polys):
                rotated = rotate_polygon(base, rotor_overlay.pivot, rotor_angle)
                patch.set_xy(rotated)
                artists.append(patch)
            for patch, base in zip(magnet_patches, magnet_base_polys):
                rotated = rotate_polygon(base, rotor_overlay.pivot, rotor_angle)
                patch.set_xy(rotated)
                artists.append(patch)
        return artists

    update(0)

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(datasets),
        blit=False,
        interval=1000 / fps,
    )

    save_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = save_path.suffix.lower()
    if suffix == ".gif":
        writer = animation.PillowWriter(fps=fps)
    else:
        writer = animation.FFMpegWriter(fps=fps)
    anim.save(str(save_path), writer=writer)

    if still_path is not None:
        still_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(still_path, dpi=150)

    if html_path is not None:
        html_path.parent.mkdir(parents=True, exist_ok=True)
        video_html = anim.to_jshtml(fps=fps)
        phase_info = """
<details>
<summary>How to read this animation</summary>
<p>The top panel shows the spatial magnetic flux density magnitude with field vectors and the live rotor overlay.
The middle gauge indicates the dominant bore field direction, while the bottom plot reports the available supply currents.</p>
</details>
"""
        html_path.write_text(
            "<html><head><meta charset='utf-8'><title>Motor field animation</title></head><body>"
            + "<h1>Motor field animation</h1>"
            + phase_info
            + video_html
            + "</body></html>",
            encoding="utf-8",
        )

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pvd", type=Path, required=True, help="VTK time-series PVD file")
    parser.add_argument("--scenario", type=Path, required=True, help="Scenario JSON used for the solve")
    parser.add_argument("--save", type=Path, required=True, help="Output animation path (e.g. MP4)")
    parser.add_argument("--fps", type=int, default=24, help="Animation frames per second (default: 24)")
    parser.add_argument("--width", type=int, default=640, help="Figure width in pixels (default: 640)")
    parser.add_argument("--html", type=Path, help="Optional HTML output embedding an interactive player")
    parser.add_argument("--frame-png", type=Path, help="Optional path for saving the first frame as a PNG image")
    parser.add_argument("--log-scale", action="store_true", help="Render |B| with a logarithmic color scale")
    parser.add_argument("--rotor", help="Rotor name to overlay (defaults to the first rotor in the scenario)")
    parser.add_argument(
        "--mechanical",
        type=Path,
        help="Mechanical trace CSV used to recover rotor angles for the overlay",
    )
    parser.add_argument(
        "--mechanical-rotor",
        help="Rotor identifier used inside the mechanical trace (defaults to --rotor)",
    )
    parser.add_argument(
        "--circuit-trace",
        type=Path,
        help="Circuit trace CSV for plotting currents when the scenario timeline lacks phase currents",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scenario = load_scenario(args.scenario)
    datasets = load_pvd(args.pvd)
    if not datasets:
        raise RuntimeError("No datasets found in PVD file")
    bore_polygon = extract_bore_polygon(scenario)
    bore_times, bore_bx, bore_by = compute_bore_series(datasets, bore_polygon)
    if len(datasets) != bore_times.size:
        raise RuntimeError("VTK frame count mismatch")

    rotor_name = args.rotor
    rotor_overlay = extract_rotor_overlay(scenario, rotor_name)
    if rotor_overlay is not None and rotor_name is None:
        rotor_name = rotor_overlay.name

    explicit_mechanical = args.mechanical is not None
    mechanical_path: Optional[Path] = args.mechanical
    if mechanical_path is None:
        mechanical_path = infer_output_path(scenario, "mechanical_trace", args.scenario)
    mechanical_series: Optional[MechanicalSeries] = None
    mechanical_key = args.mechanical_rotor or rotor_name
    if mechanical_path is not None:
        resolved_mechanical = resolve_optional_path(args.scenario.parent, str(mechanical_path))
        if resolved_mechanical.exists():
            if mechanical_key is None:
                raise RuntimeError("Mechanical trace provided but rotor name is unknown")
            mechanical_series = load_mechanical_trace(resolved_mechanical, mechanical_key)
        elif explicit_mechanical:
            raise RuntimeError(f"Mechanical trace '{resolved_mechanical}' not found")
    if mechanical_series is None and rotor_name is not None:
        mechanical_series = extract_timeline_mechanical(scenario, rotor_name)

    rotor_angles = sample_rotor_angles(mechanical_series, bore_times)

    circuit_path: Optional[Path] = args.circuit_trace
    if circuit_path is not None:
        resolved_circuit = resolve_optional_path(args.scenario.parent, str(circuit_path))
        if not resolved_circuit.exists():
            raise RuntimeError(f"Circuit trace '{resolved_circuit}' not found")
        circuit_path = resolved_circuit
    current_series = load_current_series(scenario, circuit_path)

    build_animation(
        args.save,
        args.html,
        args.frame_png,
        scenario,
        datasets,
        current_series,
        bore_times,
        bore_bx,
        bore_by,
        bore_polygon,
        rotor_overlay,
        rotor_angles,
        fps=args.fps,
        width=args.width,
        log_scale=args.log_scale,
    )


if __name__ == "__main__":
    main()
