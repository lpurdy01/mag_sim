"""Create animations for the three-phase stator rotating field."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import vtk  # type: ignore
from matplotlib import patheffects
from matplotlib.colors import LogNorm
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
    log_scale: bool,
) -> None:
    x_coords, y_coords, field_frames = load_field_frames(datasets)
    bore_magnitude = np.hypot(bore_bx, bore_by)
    max_mag = float(np.max(bore_magnitude)) if bore_magnitude.size > 0 else 1.0
    bore_angles = np.arctan2(bore_by, bore_bx)
    outlines = extract_outline_polygons(scenario)
    slot_annotations = extract_slot_annotations(scenario)

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

    quiver_step = max(1, int(max(first_mag.shape) / 30))
    quiver = ax_field.quiver(
        x_coords[::quiver_step, ::quiver_step],
        y_coords[::quiver_step, ::quiver_step],
        field_frames[0]["bx"][::quiver_step, ::quiver_step],
        field_frames[0]["by"][::quiver_step, ::quiver_step],
        color="white",
        scale=None,
        width=0.005,
    )

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

    phase_colors = {"A": "tab:red", "B": "tab:green", "C": "tab:blue"}
    for slot in slot_annotations:
        vertices = slot.get("vertices")
        if not isinstance(vertices, np.ndarray) or vertices.size == 0:
            continue
        cx, cy = polygon_centroid(vertices)
        text = ax_field.text(
            cx,
            cy,
            str(slot.get("label", "")),
            color=phase_colors.get(slot.get("phase", ""), "white"),
            fontsize=9,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=6,
        )
        text.set_path_effects([patheffects.withStroke(linewidth=3.0, foreground="black")])

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

    ax_currents.plot(phase_times, currents_a, label="Ia", color="tab:red")
    ax_currents.plot(phase_times, currents_b, label="Ib", color="tab:green")
    ax_currents.plot(phase_times, currents_c, label="Ic", color="tab:blue")
    ax_currents.set_xlabel("Time [s]")
    ax_currents.set_ylabel("Current [A]")
    ax_currents.legend(loc="upper right")
    ax_currents.grid(True, alpha=0.3)
    marker_line = ax_currents.axvline(phase_times[0] if phase_times.size else 0.0, color="k", linestyle="--")

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
        return [im, quiver, arrow, angle_text, marker_line]

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
<p>The top panel shows the spatial magnetic flux density magnitude with field vectors.
The middle gauge indicates the dominant bore field direction, while the bottom plot overlays the three balanced phase currents.</p>
</details>
"""
        html_path.write_text(
            "<html><head><meta charset='utf-8'><title>Three-phase stator field animation</title></head><body>"
            + "<h1>Three-phase stator rotating field</h1>"
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
        args.html,
        args.frame_png,
        scenario,
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
        log_scale=args.log_scale,
    )


if __name__ == "__main__":
    main()
