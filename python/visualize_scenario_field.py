#!/usr/bin/env python3
"""Visualise a scenario field map and overlay wire positions."""

from __future__ import annotations

import argparse
import json
import pathlib
from dataclasses import dataclass
from typing import List, Tuple

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
    )


if __name__ == "__main__":
    main()
