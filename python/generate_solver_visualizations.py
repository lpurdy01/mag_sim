#!/usr/bin/env python3
"""Generate warm-start and progress snapshot animations for the solver stack."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple
from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter


@dataclass
class HistorySeries:
    iterations: np.ndarray
    residuals: np.ndarray
    elapsed: np.ndarray


@dataclass
class SnapshotRecord:
    iteration: int
    x: np.ndarray
    y: np.ndarray
    field: np.ndarray


def run_capture(binary: Path) -> None:
    if not binary.exists():
        raise FileNotFoundError(f"Capture tool not found: {binary}")
    print(f"[capture] Running {binary} …", flush=True)
    subprocess.run([str(binary)], check=True)


def load_history_csv(path: Path) -> HistorySeries:
    iterations: List[int] = []
    residuals: List[float] = []
    elapsed: List[float] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            iterations.append(int(float(row["iter"])))
            residuals.append(float(row["rel_residual"]))
            elapsed.append(float(row["elapsed_seconds"]))
    if not iterations:
        raise ValueError(f"History CSV {path} is empty")
    return HistorySeries(np.asarray(iterations), np.asarray(residuals), np.asarray(elapsed))


def load_snapshot_csv(path: Path) -> SnapshotRecord:
    xs: List[float] = []
    ys: List[float] = []
    values: List[float] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))
            values.append(float(row["Az"]))
    if not values:
        raise ValueError(f"Snapshot CSV {path} is empty")
    unique_x = np.unique(np.asarray(xs))
    unique_y = np.unique(np.asarray(ys))
    grid = np.asarray(values).reshape(unique_y.size, unique_x.size)
    iteration = extract_iteration_from_name(path.stem)
    return SnapshotRecord(iteration=iteration, x=unique_x, y=unique_y, field=grid)


@dataclass
class ScenarioOverlay:
    origin_x: float
    origin_y: float
    ring_polygons: List[np.ndarray]
    magnets: List[np.ndarray]
    wires: List[Tuple[float, float, float]]  # (x, y, radius)


def load_scenario_overlay(path: Path) -> ScenarioOverlay:
    with path.open() as handle:
        data = json.load(handle)

    domain = data["domain"]
    origin_x = -0.5 * float(domain["Lx"])
    origin_y = -0.5 * float(domain["Ly"])

    ring_polygons: List[np.ndarray] = []
    for region in data.get("regions", []):
        if region.get("type") == "polygon" and region.get("material") == "iron":
            vertices = np.asarray(region.get("vertices", []), dtype=float)
            if vertices.size:
                ring_polygons.append(np.vstack([vertices, vertices[0]]))

    magnets: List[np.ndarray] = []
    for magnet in data.get("magnet_regions", []):
        vertices = np.asarray(magnet.get("vertices", []), dtype=float)
        if vertices.size:
            magnets.append(np.vstack([vertices, vertices[0]]))

    wires: List[Tuple[float, float, float]] = []
    for source in data.get("sources", []):
        if source.get("type") == "wire":
            wires.append((float(source.get("x", 0.0)),
                          float(source.get("y", 0.0)),
                          float(source.get("radius", 0.0))))

    return ScenarioOverlay(origin_x=origin_x,
                           origin_y=origin_y,
                           ring_polygons=ring_polygons,
                           magnets=magnets,
                           wires=wires)
def extract_iteration_from_name(stem: str) -> int:
    match = re.search(r"(\d+)$", stem)
    if not match:
        raise ValueError(f"Unable to infer iteration from filename '{stem}'")
    return int(match.group(1))


def create_warm_start_animation(cold: HistorySeries,
                                warm: HistorySeries,
                                frame_index: int,
                                output_path: Path,
                                image_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(f"Frame {frame_index}: Warm-Start vs Cold-Start Convergence Comparison", fontsize=12)
    ax.set_xlabel("Iteration number", fontsize=11)
    ax.set_ylabel("Relative residual (dimensionless)", fontsize=11)
    ax.set_yscale("log")
    max_iter = max(cold.iterations.max(), warm.iterations.max())
    ax.set_xlim(0, max_iter * 1.02)
    min_residual = min(cold.residuals.min(), warm.residuals.min())
    max_residual = max(cold.residuals.max(), warm.residuals.max())
    ax.set_ylim(min_residual * 0.9, max_residual * 1.1)
    ax.grid(True, alpha=0.3, linestyle=':')

    cold_line, = ax.plot([], [], label="Cold start (no initial guess)", color="#1f77b4", linewidth=2)
    warm_line, = ax.plot([], [], label="Warm start (prior solution)", color="#ff7f0e", linewidth=2)
    ax.legend(fontsize=10, loc='upper right', framealpha=0.9)

    total_frames = max(len(cold.iterations), len(warm.iterations))

    def frame_data(series: HistorySeries, idx: int) -> Tuple[np.ndarray, np.ndarray]:
        limit = min(idx + 1, len(series.iterations))
        return series.iterations[:limit], series.residuals[:limit]

    def update(frame: int):
        cold_x, cold_y = frame_data(cold, frame)
        warm_x, warm_y = frame_data(warm, frame)
        cold_line.set_data(cold_x, cold_y)
        warm_line.set_data(warm_x, warm_y)
        return cold_line, warm_line

    animation = FuncAnimation(fig, update, frames=total_frames, interval=60, blit=False)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    animation.save(output_path, writer=PillowWriter(fps=20))

    # Save a static comparison plot as well.
    cold_line.set_data(cold.iterations, cold.residuals)
    warm_line.set_data(warm.iterations, warm.residuals)
    fig.savefig(image_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def create_snapshot_animation(history: HistorySeries,
                              snapshots: Sequence[SnapshotRecord],
                              output_path: Path,
                              overlay: ScenarioOverlay,
                              solver_label: str) -> None:
    if not snapshots:
        raise ValueError("No snapshots available for animation")

    sorted_snaps = sorted(snapshots, key=lambda snap: snap.iteration)
    global_min = min(np.min(snap.field) for snap in sorted_snaps)
    global_max = max(np.max(snap.field) for snap in sorted_snaps)

    fig, (ax_field, ax_curve) = plt.subplots(1, 2, figsize=(12.5, 5))
    fig.subplots_adjust(wspace=0.38)

    initial = sorted_snaps[0]
    extent = [initial.x.min() + overlay.origin_x,
              initial.x.max() + overlay.origin_x,
              initial.y.min() + overlay.origin_y,
              initial.y.max() + overlay.origin_y]

    im = ax_field.imshow(initial.field,
                         extent=extent,
                         origin="lower",
                         cmap="viridis",
                         vmin=global_min,
                         vmax=global_max,
                         aspect="equal")
    ax_field.set_title(f"Magnetic Vector Potential at Iteration {initial.iteration}", fontsize=11)
    ax_field.set_xlabel("x position (m)", fontsize=10)
    ax_field.set_ylabel("y position (m)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax_field, shrink=0.85, label="Az (Wb/m)")
    cbar.ax.tick_params(labelsize=9)
    
    # Add geometry outlines derived from the scenario definition
    if overlay.ring_polygons:
        ax_field.plot(overlay.ring_polygons[0][:, 0],
                      overlay.ring_polygons[0][:, 1],
                      'r--', linewidth=1.5, alpha=0.8, label='Iron ring')
        for poly in overlay.ring_polygons[1:]:
            ax_field.plot(poly[:, 0], poly[:, 1], 'r--', linewidth=1.5, alpha=0.8)

    for idx, (wx, wy, radius) in enumerate(overlay.wires):
        circle = plt.Circle((wx, wy), radius, color='orange', fill=False,
                             linewidth=1.2, alpha=0.85, linestyle='-')
        if idx == 0:
            circle.set_label('Excitation wires')
        ax_field.add_patch(circle)

    if overlay.magnets:
        ax_field.plot(overlay.magnets[0][:, 0], overlay.magnets[0][:, 1],
                      'b-', linewidth=2, alpha=0.85, label='Permanent magnet')
        for magnet in overlay.magnets[1:]:
            ax_field.plot(magnet[:, 0], magnet[:, 1], 'b-', linewidth=2, alpha=0.85)
    
    ax_field.legend(loc='upper right', fontsize=8, framealpha=0.9)

    ax_curve.set_title(f"{solver_label} Solver Convergence History", fontsize=11)
    ax_curve.set_xlabel("Iteration number", fontsize=10)
    ax_curve.set_ylabel("Relative residual (dimensionless)", fontsize=10)
    ax_curve.set_yscale("log")
    ax_curve.set_xlim(0, history.iterations.max() * 1.02)
    ax_curve.set_ylim(history.residuals.min() * 0.9, history.residuals.max() * 1.1)
    ax_curve.grid(True, alpha=0.3, linestyle=':')

    progress_line, = ax_curve.plot([], [], color="#1f77b4")
    marker, = ax_curve.plot([], [], marker="o", color="#d62728")

    def update(frame: int):
        snapshot = sorted_snaps[frame]
        im.set_data(snapshot.field)
        im.set_extent([snapshot.x.min() + overlay.origin_x,
                       snapshot.x.max() + overlay.origin_x,
                       snapshot.y.min() + overlay.origin_y,
                       snapshot.y.max() + overlay.origin_y])
        ax_field.set_title(f"Magnetic Vector Potential (Iteration {snapshot.iteration})")

        mask = history.iterations <= snapshot.iteration
        progress_line.set_data(history.iterations[mask], history.residuals[mask])
        if np.any(mask):
            current_iter = history.iterations[mask][-1]
            current_res = history.residuals[mask][-1]
            marker.set_data([current_iter], [current_res])
        else:
            marker.set_data([], [])
        return im, progress_line, marker

    animation = FuncAnimation(fig, update, frames=len(sorted_snaps), interval=400, blit=False)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    animation.save(output_path, writer=PillowWriter(fps=2))
    plt.close(fig)


def create_convergence_overlay_plot(histories: Dict[str, HistorySeries],
                                    output_path: Path,
                                    iter_limit: int = 2000) -> None:
    plt.figure(figsize=(7.5, 5))
    for label, series in histories.items():
        mask = series.iterations <= iter_limit
        if not np.any(mask):
            mask = series.iterations <= series.iterations.max()
        plt.plot(series.iterations[mask], series.residuals[mask], label=label, linewidth=2)
    plt.yscale("log")
    plt.xlabel("Iteration number", fontsize=11)
    plt.ylabel("Relative residual (dimensionless)", fontsize=11)
    plt.title("CG vs SOR Convergence (first 2000 iterations)", fontsize=12)
    plt.grid(True, linestyle=":", alpha=0.3)
    plt.xlim(0, iter_limit)
    plt.legend(fontsize=10)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close()


def discover_frame_history(base_dir: Path, frame_index: int, warm: bool) -> Path:
    subdir = "warm" if warm else "cold"
    path = base_dir / "warm_start" / subdir / f"frame_{frame_index}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing history for frame {frame_index} ({'warm' if warm else 'cold'}): {path}")
    return path


def discover_snapshot_files(base_dir: Path, solver_label: str) -> Tuple[HistorySeries, List[SnapshotRecord]]:
    history_path = base_dir / "progress_snapshots" / solver_label / "history.csv"
    if not history_path.exists():
        raise FileNotFoundError(f"Snapshot history not found: {history_path}")
    history = load_history_csv(history_path)

    snapshot_dir = base_dir / "progress_snapshots" / solver_label / "snapshots"
    snapshot_paths = sorted(snapshot_dir.glob("snapshot_iter_*.csv"))
    snapshots = [load_snapshot_csv(path) for path in snapshot_paths]
    return history, snapshots


def run_solver_benchmark_series(build_dir: Path,
                                output_dir: Path,
                                grid_sizes: Iterable[int] | Mapping[str, Iterable[int]],
                                repeats: int = 3,
                                max_iters: int = 20000,
                                tol: float = 1e-6) -> List[Dict[str, float]]:
    benchmark_exe = build_dir / "solver_benchmark"
    if not benchmark_exe.exists():
        raise FileNotFoundError(f"solver_benchmark not found: {benchmark_exe}")

    output_dir.mkdir(parents=True, exist_ok=True)
    results: List[Dict[str, float]] = []

    if isinstance(grid_sizes, Mapping):
        solver_size_map = {key: list(value) for key, value in grid_sizes.items()}
    else:
        shared = list(grid_sizes)
        solver_size_map = {"cg": shared, "sor": shared}

    for solver_label in ("cg", "sor"):
        sizes = solver_size_map.get(solver_label)
        if not sizes:
            continue
        csv_path = output_dir / f"{solver_label}_benchmark.csv"
        if csv_path.exists():
            csv_path.unlink()
        for size in sizes:
            cmd = [
                str(benchmark_exe),
                "--nx", str(size),
                "--ny", str(size),
                "--max-iters", str(max_iters),
                "--tol", str(tol),
                "--repeats", str(repeats),
                "--solver", solver_label,
                "--csv", str(csv_path),
            ]
            if solver_label == "sor":
                cmd.extend(["--omega", "1.8"])
            subprocess.run(cmd, check=True)

        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                entry = {
                    "solver": solver_label.upper(),
                    "nx": float(row["nx"]),
                    "ny": float(row["ny"]),
                    "cells": float(row["cells"]),
                    "iters": float(row["iters"]),
                    "avg_ms": float(row["avg_ms"]),
                    "updates_per_second": float(row["updates_per_second"]),
                }
                results.append(entry)

    return results


def create_scaling_plot_and_report(results: List[Dict[str, float]],
                                   output_dir: Path,
                                   frames_per_cycle: int = 120) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    solvers = sorted({entry["solver"] for entry in results})
    size_values = sorted({int(entry["nx"]) for entry in results})

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].set_title("Solve Time Scaling vs Grid Size", fontsize=12)
    axes[0].set_xlabel("Grid edge cells (N)", fontsize=11)
    axes[0].set_ylabel("Average solve time (s)", fontsize=11)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")

    summary_lines: List[str] = []

    scaling_stats: Dict[str, Tuple[float, float]] = {}
    for solver in solvers:
        solver_entries = [entry for entry in results if entry["solver"] == solver]
        solver_entries.sort(key=lambda e: e["nx"])
        ns = np.array([entry["nx"] for entry in solver_entries])
        avg_seconds = np.array([entry["avg_ms"] for entry in solver_entries]) / 1000.0
        axes[0].plot(ns, avg_seconds, marker="o", linewidth=2, label=f"{solver}")

        log_n = np.log(ns)
        log_t = np.log(avg_seconds)
        slope, intercept = np.polyfit(log_n, log_t, 1)
        scaling_stats[solver] = (slope, intercept)
        summary_lines.append(f"{solver} scaling exponent ~ N^{slope:.2f}")

        fit_t = np.exp(intercept) * ns**slope
        axes[0].plot(ns, fit_t, linestyle="--", linewidth=1.2, alpha=0.6)

    axes[0].grid(True, which="both", linestyle=":", alpha=0.3)
    axes[0].legend(fontsize=10)

    # Per-cycle timing (assuming frames_per_cycle solves)
    axes[1].set_title(f"CPU Time per {frames_per_cycle}-frame Cycle", fontsize=12)
    axes[1].set_xlabel("Grid edge cells (N)", fontsize=11)
    axes[1].set_ylabel("Seconds per cycle", fontsize=11)
    axes[1].set_xscale("log")
    axes[1].grid(True, which="both", linestyle=":", alpha=0.3)

    for solver in solvers:
        solver_entries = [entry for entry in results if entry["solver"] == solver]
        solver_entries.sort(key=lambda e: e["nx"])
        ns = np.array([entry["nx"] for entry in solver_entries])
        per_cycle_seconds = np.array([entry["avg_ms"] for entry in solver_entries]) / 1000.0 * frames_per_cycle
        axes[1].plot(ns, per_cycle_seconds, marker="o", linewidth=2, label=solver)

    axes[1].legend(fontsize=10)

    plot_path = output_dir / "solver_scaling.png"
    fig.tight_layout()
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)

    # Detailed text summary
    summary_path = output_dir / "solver_scaling_summary.txt"
    header = "GridN, Cells, Solver, AvgSolveMs, AvgSolvePerCycleSec, Iterations, UpdatesPerSecond\n"
    lines = [header]
    for size in size_values:
        for solver in solvers:
            entry = next((e for e in results if e["solver"] == solver and int(e["nx"]) == size), None)
            if entry is None:
                continue
            avg_sec = entry["avg_ms"] / 1000.0
            per_cycle_sec = avg_sec * frames_per_cycle
            lines.append(
                f"{size},{int(entry['cells'])},{solver},{entry['avg_ms']:.3f},{per_cycle_sec:.3f},"
                f"{entry['iters']:.0f},{entry['updates_per_second']:.3e}\n"
            )

    lines.append("\nScaling fits (solve time ≈ a * N^k):\n")
    for solver, (slope, intercept) in scaling_stats.items():
        coeff = np.exp(intercept)
        lines.append(f"  {solver}: k={slope:.2f}, a={coeff:.3e}\n")

    summary_lines_text = "\n".join(summary_lines)
    lines.append("\n" + summary_lines_text + "\n")

    with summary_path.open("w") as handle:
        handle.writelines(lines)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build"), help="CMake build directory")
    parser.add_argument("--skip-capture", action="store_true", help="Reuse existing CSVs instead of regenerating")
    parser.add_argument("--frame", type=int, default=2, help="Timeline frame index to visualise for warm-start animation")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs") / "visualization" / "animations",
                        help="Directory for generated animations")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    build_dir: Path = args.build_dir
    capture_binary = build_dir / "solver_visualization_capture"

    if not args.skip_capture:
        run_capture(capture_binary)
    else:
        print("[capture] Skipping dataset regeneration", flush=True)

    dataset_root = Path("outputs") / "visualization"
    scenario_path = Path("inputs") / "iron_ring_magnet_viz.json"
    cold_path = discover_frame_history(dataset_root, args.frame, warm=False)
    warm_path = discover_frame_history(dataset_root, args.frame, warm=True)

    cold_history = load_history_csv(cold_path)
    warm_history = load_history_csv(warm_path)

    animations_dir = args.output_dir
    warm_animation_path = animations_dir / f"warm_start_frame_{args.frame}.gif"
    warm_image_path = animations_dir / f"warm_start_frame_{args.frame}.png"
    print(f"[viz] Creating warm-start animation at {warm_animation_path}")
    create_warm_start_animation(cold_history, warm_history, args.frame, warm_animation_path, warm_image_path)

    overlay = load_scenario_overlay(scenario_path)
    solver_variants = [
        ("cg", "CG"),
        ("sor", "SOR"),
    ]

    overlay_histories: Dict[str, HistorySeries] = {}

    for label, pretty in solver_variants:
        try:
            history, snapshots = discover_snapshot_files(dataset_root, label)
        except FileNotFoundError as exc:
            print(f"[viz] Skipping {pretty} snapshot animation: {exc}")
            continue
        snapshot_animation_path = animations_dir / f"progress_snapshots_{label}.gif"
        print(f"[viz] Creating progress snapshot animation at {snapshot_animation_path}")
        create_snapshot_animation(history, snapshots, snapshot_animation_path, overlay, pretty)
        overlay_histories[pretty] = history

    if overlay_histories:
        overlay_plot_path = animations_dir / "convergence_overlay.png"
        create_convergence_overlay_plot(overlay_histories, overlay_plot_path)

    benchmark_dir = Path("outputs") / "visualization" / "benchmark"
    solver_size_map = {
        "cg": [65, 97, 129, 193, 257, 385, 513, 769, 961],
        "sor": [65, 97, 129],
    }
    benchmark_results = run_solver_benchmark_series(build_dir, benchmark_dir, solver_size_map)
    create_scaling_plot_and_report(benchmark_results, benchmark_dir)

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

