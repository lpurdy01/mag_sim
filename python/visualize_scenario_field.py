#!/usr/bin/env python3
"""Visualise a scenario field map and overlay wire positions."""

from __future__ import annotations

import argparse
import json
import pathlib
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - optional dependency guard
    raise SystemExit(
        "numpy is required for visualization. Install it via 'pip install numpy'."
    ) from exc


@dataclass
class Wire:
    x: float
    y: float
    radius: float
    current: float


MU0 = 4.0e-7 * np.pi


def compute_domain_bounds(spec: dict) -> Tuple[float, float, float, float]:
    domain = spec.get("domain", {})
    lx = float(domain.get("Lx", 0.0))
    ly = float(domain.get("Ly", 0.0))
    if lx <= 0.0 or ly <= 0.0:
        raise ValueError("Scenario domain must define positive Lx and Ly to draw boundaries")
    xmin = -0.5 * lx
    xmax = 0.5 * lx
    ymin = -0.5 * ly
    ymax = 0.5 * ly
    return xmin, xmax, ymin, ymax


def line_rectangle_intersections(nx: float, ny: float, offset: float,
                                 bounds: Tuple[float, float, float, float]) -> List[Tuple[float, float]]:
    xmin, xmax, ymin, ymax = bounds
    points: List[Tuple[float, float]] = []

    def add_point(x: float, y: float) -> None:
        for px, py in points:
            if abs(px - x) < 1e-9 and abs(py - y) < 1e-9:
                return
        points.append((x, y))

    if abs(ny) > 1e-12:
        for x in (xmin, xmax):
            y = -(nx * x + offset) / ny
            if ymin - 1e-9 <= y <= ymax + 1e-9:
                add_point(x, y)
    if abs(nx) > 1e-12:
        for y in (ymin, ymax):
            x = -(ny * y + offset) / nx
            if xmin - 1e-9 <= x <= xmax + 1e-9:
                add_point(x, y)
    return points[:2]


def gather_boundaries(spec: dict) -> Dict[str, List]:
    bounds = compute_domain_bounds(spec)
    regions = spec.get("regions", [])
    materials = {entry.get("name"): float(entry.get("mu_r", 1.0)) for entry in spec.get("materials", [])}

    halfspaces: List[Tuple[Tuple[float, float], Tuple[float, float], float]] = []
    polygons: List[Tuple[List[float], List[float], str]] = []
    for region in regions:
        rtype = region.get("type")
        if rtype == "halfspace":
            normal = region.get("normal", [0.0, 0.0])
            if len(normal) != 2:
                continue
            nx, ny = float(normal[0]), float(normal[1])
            offset = float(region.get("offset", 0.0))
            pts = line_rectangle_intersections(nx, ny, offset, bounds)
            if len(pts) == 2:
                material = str(region.get("material", ""))
                mu = materials.get(material, float("nan"))
                halfspaces.append((pts[0], pts[1], mu))
        elif rtype == "polygon":
            vertices = region.get("vertices", [])
            xs = [float(v[0]) for v in vertices if isinstance(v, list) and len(v) == 2]
            ys = [float(v[1]) for v in vertices if isinstance(v, list) and len(v) == 2]
            if len(xs) >= 3:
                material = str(region.get("material", ""))
                polygons.append((xs, ys, material))

    magnet_polys: List[Tuple[List[float], List[float]]] = []
    magnet_rects: List[Tuple[float, float, float, float]] = []
    for region in spec.get("magnet_regions", []):
        rtype = region.get("type")
        if rtype == "polygon":
            vertices = region.get("vertices", [])
            xs = [float(v[0]) for v in vertices if isinstance(v, list) and len(v) == 2]
            ys = [float(v[1]) for v in vertices if isinstance(v, list) and len(v) == 2]
            if len(xs) >= 3:
                magnet_polys.append((xs, ys))
        elif rtype == "rect":
            xrange = region.get("x_range", [0.0, 0.0])
            yrange = region.get("y_range", [0.0, 0.0])
            if len(xrange) == 2 and len(yrange) == 2:
                magnet_rects.append((float(xrange[0]), float(xrange[1]), float(yrange[0]), float(yrange[1])))

    return {
        "halfspaces": halfspaces,
        "polygons": polygons,
        "magnet_polygons": magnet_polys,
        "magnet_rects": magnet_rects,
    }


def compute_interface_reference_field(
    spec: dict, grid_x: np.ndarray, grid_y: np.ndarray
) -> Optional[Dict[str, np.ndarray]]:
    regions = [region for region in spec.get("regions", []) if region.get("type") == "halfspace"]
    if len(regions) != 2:
        return None

    normals = [region.get("normal", [0.0, 0.0]) for region in regions]
    offsets = [float(region.get("offset", 0.0)) for region in regions]

    try:
        nx_vals = [float(n[0]) for n in normals]
        ny_vals = [float(n[1]) for n in normals]
    except (TypeError, ValueError):
        return None

    if any(abs(ny) > 1e-6 for ny in ny_vals):
        return None

    left_index = None
    right_index = None
    for idx, nx in enumerate(nx_vals):
        if nx > 0.0:
            left_index = idx
        elif nx < 0.0:
            right_index = idx
    if left_index is None or right_index is None:
        return None

    materials = {entry.get("name"): float(entry.get("mu_r", 1.0)) for entry in spec.get("materials", [])}
    try:
        mu_left = MU0 * materials[regions[left_index]["material"]]
        mu_right = MU0 * materials[regions[right_index]["material"]]
    except KeyError:
        return None

    sources = [src for src in spec.get("sources", []) if src.get("type") == "wire"]
    if not sources:
        return None
    wire = sources[0]
    wire_x = float(wire.get("x", 0.0))
    wire_y = float(wire.get("y", 0.0))
    current = float(wire.get("I", 0.0))

    interface_x = -offsets[left_index] / nx_vals[left_index]
    image_x = 2.0 * interface_x - wire_x
    rho = (mu_left - mu_right) / (mu_left + mu_right)
    tau = (2.0 * mu_right) / (mu_left + mu_right)

    dx_real = grid_x - wire_x
    dy_real = grid_y - wire_y
    r2_real = dx_real * dx_real + dy_real * dy_real

    dx_image = grid_x - image_x
    dy_image = grid_y - wire_y
    r2_image = dx_image * dx_image + dy_image * dy_image

    with np.errstate(divide="ignore", invalid="ignore"):
        coeff_real_left = (mu_left * current) / (2.0 * np.pi * r2_real)
        coeff_image = (mu_left * rho * current) / (2.0 * np.pi * r2_image)
        coeff_transmitted = (mu_right * tau * current) / (2.0 * np.pi * r2_real)

    coeff_real_left = np.where(r2_real > 0.0, coeff_real_left, 0.0)
    coeff_image = np.where(r2_image > 0.0, coeff_image, 0.0)
    coeff_transmitted = np.where(r2_real > 0.0, coeff_transmitted, 0.0)

    boundary_eval = (
        nx_vals[left_index] * grid_x + ny_vals[left_index] * grid_y + offsets[left_index]
    )
    left_mask = boundary_eval < 0.0

    bx = np.zeros_like(grid_x)
    by = np.zeros_like(grid_y)

    bx[left_mask] = -(coeff_real_left[left_mask] * dy_real[left_mask] + coeff_image[left_mask] * dy_image[left_mask])
    by[left_mask] = coeff_real_left[left_mask] * dx_real[left_mask] + coeff_image[left_mask] * dx_image[left_mask]

    right_mask = ~left_mask
    bx[right_mask] = -coeff_transmitted[right_mask] * dy_real[right_mask]
    by[right_mask] = coeff_transmitted[right_mask] * dx_real[right_mask]

    bmag = np.hypot(bx, by)
    return {"Bx": bx, "By": by, "Bmag": bmag}
def load_scenario(path: pathlib.Path) -> Tuple[dict, List[Wire]]:
    with path.open("r", encoding="utf-8") as handle:
        spec = json.load(handle)

    wires: List[Wire] = []
    for source in spec.get("sources", []):
        if source.get("type") != "wire":
            continue
        wires.append(
            Wire(
                x=float(source.get("x", 0.0)),
                y=float(source.get("y", 0.0)),
                radius=float(source.get("radius", 0.0)),
                current=float(source.get("I", 0.0)),
            )
        )
    return spec, wires


def resolve_field_map_path(spec: dict, args: argparse.Namespace) -> pathlib.Path:
    if args.field_map:
        return args.field_map

    outputs = spec.get("outputs", [])
    field_maps = [out for out in outputs if out.get("type") == "field_map"]
    if not field_maps:
        raise SystemExit("Scenario defines no field_map outputs; specify --field-map explicitly.")

    if args.field_map_id:
        for out in field_maps:
            if out.get("id") == args.field_map_id:
                return pathlib.Path(out.get("path", f"outputs/{out['id']}.csv"))
        raise SystemExit(f"No field_map output with id '{args.field_map_id}' in scenario.")

    # Default to the first declared field map
    first = field_maps[0]
    return pathlib.Path(first.get("path", f"outputs/{first['id']}.csv"))


def load_field_csv(path: pathlib.Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xs = []
    ys = []
    bxs = []
    bys = []
    bmags = []

    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split(",")
        expected = ["x", "y", "Bx", "By", "Bmag"]
        if header != expected:
            raise SystemExit(f"Unexpected CSV header {header}; expected {expected}.")
        for line in handle:
            if not line.strip():
                continue
            parts = line.strip().split(",")
            if len(parts) != 5:
                raise SystemExit(f"Malformed row in CSV: {line.strip()!r}")
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
            bxs.append(float(parts[2]))
            bys.append(float(parts[3]))
            bmags.append(float(parts[4]))

    return (
        np.asarray(xs),
        np.asarray(ys),
        np.asarray(bxs),
        np.asarray(bys),
        np.asarray(bmags),
    )


def reshape_field(xs: np.ndarray, ys: np.ndarray, values: Tuple[np.ndarray, ...]) -> Tuple[np.ndarray, np.ndarray, Tuple[np.ndarray, ...]]:
    unique_x = np.unique(xs)
    unique_y = np.unique(ys)
    nx = unique_x.size
    ny = unique_y.size
    if xs.size != nx * ny:
        raise SystemExit("CSV data cannot be reshaped onto a regular grid")

    order = np.argsort(ys, kind="mergesort")
    xs = xs[order]
    ys = ys[order]
    reshaped = []
    for value in values:
        reshaped.append(value[order].reshape(ny, nx))
    grid_x = xs.reshape(ny, nx)
    grid_y = ys.reshape(ny, nx)
    return grid_x, grid_y, tuple(reshaped)


def plot_field(
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    bx: np.ndarray,
    by: np.ndarray,
    bmag: np.ndarray,
    wires: List[Wire],
    quiver_skip: int,
    save: pathlib.Path | None,
    title: str,
    color_scale: str,
    log_floor: float,
    draw_boundaries: bool,
    boundaries: Dict[str, List],
    show_vectors: bool,
    streamlines: bool,
    analytic: Optional[Dict[str, np.ndarray]],
    analytic_levels: int,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    display_mag = bmag
    norm = None
    if color_scale == "log":
        if log_floor <= 0.0:
            raise SystemExit("--log-floor must be positive when using log colour scale")
        positive = bmag[np.isfinite(bmag) & (bmag > 0.0)]
        if positive.size == 0:
            raise SystemExit("Cannot apply log scale because |B| has no positive samples")
        vmin = max(float(positive.min()), log_floor)
        vmax = float(positive.max())
        norm = LogNorm(vmin=vmin, vmax=vmax)
        display_mag = np.maximum(bmag, vmin)
    pcm = ax.pcolormesh(grid_x, grid_y, display_mag, shading="auto", cmap="viridis", norm=norm)
    cbar_label = "|B| [T]" if color_scale == "linear" else "|B| [T] (log scale)"
    fig.colorbar(pcm, ax=ax, label=cbar_label)

    if streamlines:
        ax.streamplot(grid_x, grid_y, bx, by, density=1.2, color="w", linewidth=0.6)

    if show_vectors:
        skip = max(1, quiver_skip)
        ax.quiver(
            grid_x[::skip, ::skip],
            grid_y[::skip, ::skip],
            bx[::skip, ::skip],
            by[::skip, ::skip],
            color="w",
            scale=0.002,
            width=0.003,
            pivot="mid",
        )

    if draw_boundaries:
        for (p0, p1, mu) in boundaries.get("halfspaces", []):
            xs = [p0[0], p1[0]]
            ys = [p0[1], p1[1]]
            ax.plot(xs, ys, color="w", linestyle="--", linewidth=1.0, alpha=0.6)
        for xs, ys, material in boundaries.get("polygons", []):
            closed_xs = xs + [xs[0]]
            closed_ys = ys + [ys[0]]
            ax.plot(closed_xs, closed_ys, color="w", linestyle="-", linewidth=1.0, alpha=0.7)
        for xs, ys in boundaries.get("magnet_polygons", []):
            closed_xs = xs + [xs[0]]
            closed_ys = ys + [ys[0]]
            ax.plot(closed_xs, closed_ys, color="tab:red", linestyle=":", linewidth=1.2, alpha=0.9)
        for xmin, xmax, ymin, ymax in boundaries.get("magnet_rects", []):
            rect_x = [xmin, xmax, xmax, xmin, xmin]
            rect_y = [ymin, ymin, ymax, ymax, ymin]
            ax.plot(rect_x, rect_y, color="tab:red", linestyle=":", linewidth=1.2, alpha=0.9)

    if analytic is not None and analytic_levels > 0:
        levels = np.linspace(np.nanmin(analytic["Bmag"]), np.nanmax(analytic["Bmag"]), analytic_levels)
        levels = levels[np.isfinite(levels)]
        if levels.size > 0:
            ax.contour(grid_x, grid_y, analytic["Bmag"], levels=levels, colors="white", linewidths=0.6, linestyles="dashed")

    for wire in wires:
        color = "tab:red" if wire.current >= 0.0 else "tab:blue"
        circle = plt.Circle((wire.x, wire.y), wire.radius, color=color, fill=False, linewidth=1.5)
        ax.add_patch(circle)
        ax.text(
            wire.x,
            wire.y + wire.radius + 0.002,
            f"I={wire.current:+.1f} A",
            color=color,
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(title)
    ax.set_xlim(grid_x.min(), grid_x.max())
    ax.set_ylim(grid_y.min(), grid_y.max())
    ax.grid(False)

    if save is not None:
        fig.savefig(save, dpi=200, bbox_inches="tight")
        print(f"Saved figure to {save}")

    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scenario",
        type=pathlib.Path,
        default=pathlib.Path("inputs/two_wire_cancel.json"),
        help="Path to the scenario JSON file.",
    )
    parser.add_argument(
        "--field-map",
        type=pathlib.Path,
        default=None,
        help="Explicit path to the field map CSV. Overrides scenario outputs if given.",
    )
    parser.add_argument(
        "--field-map-id",
        type=str,
        default=None,
        help="Name of the field_map output to use when multiple outputs are defined.",
    )
    parser.add_argument(
        "--quiver-skip",
        type=int,
        default=4,
        help="Plot every Nth vector in the quiver plot for readability.",
    )
    parser.add_argument(
        "--color-scale",
        choices=("linear", "log"),
        default="linear",
        help="Colour mapping for |B|. Use 'log' to enhance dynamic range.",
    )
    parser.add_argument(
        "--log-floor",
        type=float,
        default=1e-7,
        help="Minimum |B| value when --color-scale=log to avoid taking log of zero.",
    )
    parser.add_argument(
        "--draw-boundaries",
        action="store_true",
        help="Overlay material and magnet region boundaries on the plot.",
    )
    parser.add_argument(
        "--hide-vectors",
        action="store_true",
        help="Disable the quiver vector overlay.",
    )
    parser.add_argument(
        "--streamlines",
        action="store_true",
        help="Draw field streamlines instead of (or in addition to) vectors.",
    )
    parser.add_argument(
        "--overlay-analytic-interface",
        action="store_true",
        help=(
            "Overlay contour lines from the analytic method-of-images solution for scenarios "
            "with a single wire and planar permeability interface."
        ),
    )
    parser.add_argument(
        "--analytic-contours",
        type=int,
        default=8,
        help="Number of contour levels when plotting the analytic overlay.",
    )
    parser.add_argument(
        "--save",
        type=pathlib.Path,
        default=None,
        help="Optional path to save the rendered figure.",
    )
    args = parser.parse_args()

    if not args.scenario.exists():
        raise SystemExit(f"Scenario file not found: {args.scenario}")

    spec, wires = load_scenario(args.scenario)
    field_map_path = resolve_field_map_path(spec, args)
    if not field_map_path.exists():
        raise SystemExit(
            f"Field map CSV not found: {field_map_path}. Run the solver with --solve to generate it."
        )

    xs, ys, bxs, bys, bmags = load_field_csv(field_map_path)
    grid_x, grid_y, (bx, by, bmag) = reshape_field(xs, ys, (bxs, bys, bmags))

    boundaries: Dict[str, List] = {}
    if args.draw_boundaries or args.overlay_analytic_interface:
        try:
            boundaries = gather_boundaries(spec)
        except ValueError as exc:
            print(f"Warning: {exc}")
            boundaries = {}

    analytic = None
    if args.overlay_analytic_interface:
        analytic = compute_interface_reference_field(spec, grid_x, grid_y)
        if analytic is None:
            print("Warning: analytic interface overlay unavailable for this scenario")

    title = f"Field map: {args.scenario.name}"
    plot_field(
        grid_x,
        grid_y,
        bx,
        by,
        bmag,
        wires,
        args.quiver_skip,
        args.save,
        title,
        args.color_scale,
        args.log_floor,
        args.draw_boundaries,
        boundaries,
        not args.hide_vectors,
        args.streamlines,
        analytic,
        args.analytic_contours,
    )


if __name__ == "__main__":
    main()
