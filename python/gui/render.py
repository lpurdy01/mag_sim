"""Helpers for rendering geometry previews and field visualisations."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt  # noqa: E402  - backend configured above
from matplotlib.patches import Rectangle  # noqa: E402

from python.visualize_scenario_field import (  # noqa: E402
    Wire,
    compute_domain_bounds,
    gather_boundaries,
    load_field_csv,
    load_scenario,
    plot_field,
    reshape_field,
)


DEFAULT_LOG_FLOOR = 1e-7


def render_geometry_preview(scenario_path: Path, output_path: Path) -> Path:
    """Render a static geometry preview for the supplied scenario."""

    spec, wires = load_scenario(scenario_path)
    boundaries = {}
    try:
        boundaries = gather_boundaries(spec)
        xmin, xmax, ymin, ymax = compute_domain_bounds(spec)
    except ValueError as exc:
        raise RuntimeError(str(exc)) from exc

    fig, ax = plt.subplots(figsize=(6.0, 5.0))

    width = xmax - xmin
    height = ymax - ymin
    rect = Rectangle((xmin, ymin), width, height, fill=False, linewidth=1.6)
    ax.add_patch(rect)

    for xs, ys, _material in boundaries.get("polygons", []):
        if len(xs) < 2:
            continue
        closed_xs = xs + [xs[0]]
        closed_ys = ys + [ys[0]]
        ax.plot(closed_xs, closed_ys, color="tab:blue", linestyle="-", linewidth=1.0, alpha=0.7)

    for xs, ys in boundaries.get("magnet_polygons", []):
        if len(xs) < 2:
            continue
        closed_xs = xs + [xs[0]]
        closed_ys = ys + [ys[0]]
        ax.plot(closed_xs, closed_ys, color="tab:red", linestyle=":", linewidth=1.2, alpha=0.9)

    for xmin_m, xmax_m, ymin_m, ymax_m in boundaries.get("magnet_rects", []):
        rect_x = [xmin_m, xmax_m, xmax_m, xmin_m, xmin_m]
        rect_y = [ymin_m, ymin_m, ymax_m, ymax_m, ymin_m]
        ax.plot(rect_x, rect_y, color="tab:red", linestyle=":", linewidth=1.1, alpha=0.9)

    for wire in wires:
        color = "tab:red" if wire.current >= 0 else "tab:blue"
        circle = plt.Circle((wire.x, wire.y), wire.radius, color=color, fill=False, linewidth=1.3)
        ax.add_patch(circle)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(spec.get("name", scenario_path.name))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output_path


def render_field_map_image(
    scenario_path: Path,
    field_map_path: Path,
    *,
    output_path: Optional[Path] = None,
    vector_mode: str = "linear",
    quiver_skip: int = 4,
    color_scale: str = "linear",
    draw_boundaries: bool = True,
    streamlines: bool = False,
    log_floor: float = DEFAULT_LOG_FLOOR,
    vector_log_floor: float = DEFAULT_LOG_FLOOR,
) -> bytes | Path:
    """Render the field map with the supplied styling options."""

    spec, wires = load_scenario(scenario_path)
    xs, ys, bxs, bys, bmags = load_field_csv(field_map_path)
    grid_x, grid_y, (bx, by, bmag) = reshape_field(xs, ys, (bxs, bys, bmags))

    boundaries: Dict[str, List] = {}
    if draw_boundaries:
        try:
            boundaries = gather_boundaries(spec)
        except ValueError:
            boundaries = {}

    buffer = io.BytesIO()
    plot_field(
        grid_x,
        grid_y,
        bx,
        by,
        bmag,
        wires,
        quiver_skip,
        buffer,
        f"Field map: {scenario_path.name}",
        color_scale,
        log_floor,
        draw_boundaries,
        boundaries,
        None,
        vector_mode,
        vector_log_floor,
        streamlines,
        None,
        0,
        show=False,
    )

    data = buffer.getvalue()
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(data)
        return output_path

    return data
