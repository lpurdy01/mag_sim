"""Create animations for the three-phase stator rotating field."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import vtk  # type: ignore
from matplotlib.patches import FancyArrowPatch
from vtk.util.numpy_support import vtk_to_numpy  # type: ignore


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


def load_phase_currents(scenario: Dict[str, object]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    timeline = scenario.get("timeline", [])
    times = []
    currents_a = []
    currents_b = []
    currents_c = []
    for entry in timeline:
        if not isinstance(entry, dict):
            continue
        time = float(entry.get("t", entry.get("time", len(times))))
        phases = entry.get("phase_currents", {})
        times.append(time)
        currents_a.append(float(phases.get("A", 0.0)))
        currents_b.append(float(phases.get("B", 0.0)))
        currents_c.append(float(phases.get("C", 0.0)))
    return (
        np.asarray(times, dtype=float),
        np.asarray(currents_a, dtype=float),
        np.asarray(currents_b, dtype=float),
        np.asarray(currents_c, dtype=float),
    )


def compute_bore_series(datasets: List[Tuple[float, Path]], bore_polygon: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def build_animation(
    save_path: Path,
    datasets: List[Tuple[float, Path]],
    phase_times: np.ndarray,
    currents_a: np.ndarray,
    currents_b: np.ndarray,
    currents_c: np.ndarray,
    bore_times: np.ndarray,
    bore_bx: np.ndarray,
    bore_by: np.ndarray,
    bore_polygon: np.ndarray,
    fps: int,
    width: int,
) -> None:
    bore_magnitude = np.hypot(bore_bx, bore_by)
    max_mag = float(np.max(bore_magnitude)) if bore_magnitude.size > 0 else 1.0
    bore_angles = np.arctan2(bore_by, bore_bx)

    fig_width_in = width / 100.0
    fig_height_in = fig_width_in * 0.75
    fig, (ax_compass, ax_currents) = plt.subplots(
        2,
        1,
        figsize=(fig_width_in, fig_height_in),
        gridspec_kw={"height_ratios": [2, 1]},
    )

    bore_radius = np.max(np.linalg.norm(bore_polygon, axis=1)) if bore_polygon.size else 1.0
    circle = plt.Circle((0.0, 0.0), bore_radius, fill=False, color="tab:gray", lw=1.5)
    ax_compass.add_patch(circle)
    ax_compass.set_aspect("equal", adjustable="datalim")
    padding = bore_radius * 1.2
    ax_compass.set_xlim(-padding, padding)
    ax_compass.set_ylim(-padding, padding)
    ax_compass.axis("off")
    arrow = FancyArrowPatch((0.0, 0.0), (0.0, 0.0), color="tab:orange", arrowstyle="->", linewidth=2.0)
    ax_compass.add_patch(arrow)
    angle_text = ax_compass.text(0.02, 0.95, "", transform=ax_compass.transAxes, ha="left", va="top")

    ax_currents.plot(phase_times, currents_a, label="Ia", color="tab:red")
    ax_currents.plot(phase_times, currents_b, label="Ib", color="tab:green")
    ax_currents.plot(phase_times, currents_c, label="Ic", color="tab:blue")
    ax_currents.set_xlabel("Time [s]")
    ax_currents.set_ylabel("Current [A]")
    ax_currents.legend(loc="upper right")
    ax_currents.grid(True, alpha=0.3)
    marker_line = ax_currents.axvline(phase_times[0] if phase_times.size else 0.0, color="k", linestyle="--")

    def update(frame_idx: int) -> List[object]:
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
        return [arrow, angle_text, marker_line]

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pvd", type=Path, required=True, help="VTK time-series PVD file")
    parser.add_argument("--scenario", type=Path, required=True, help="Scenario JSON used for the solve")
    parser.add_argument("--save", type=Path, required=True, help="Output animation path (e.g. MP4)")
    parser.add_argument("--fps", type=int, default=24, help="Animation frames per second (default: 24)")
    parser.add_argument("--width", type=int, default=640, help="Figure width in pixels (default: 640)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    scenario = load_scenario(args.scenario)
    datasets = load_pvd(args.pvd)
    if not datasets:
        raise RuntimeError("No datasets found in PVD file")
    bore_polygon = extract_bore_polygon(scenario)
    phase_times, ia, ib, ic = load_phase_currents(scenario)
    bore_times, bore_bx, bore_by = compute_bore_series(datasets, bore_polygon)
    if phase_times.size != bore_times.size:
        raise RuntimeError("Timeline length and VTK frame count mismatch")
    build_animation(
        args.save,
        datasets,
        phase_times,
        ia,
        ib,
        ic,
        bore_times,
        bore_bx,
        bore_by,
        bore_polygon,
        fps=args.fps,
        width=args.width,
    )


if __name__ == "__main__":
    main()
