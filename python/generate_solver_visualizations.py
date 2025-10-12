#!/usr/bin/env python3
"""Generate warm-start and progress snapshot animations for the solver stack."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

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
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_title(f"Frame {frame_index}: warm-start vs cold-start convergence")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Relative residual")
    ax.set_yscale("log")
    max_iter = max(cold.iterations.max(), warm.iterations.max())
    ax.set_xlim(0, max_iter * 1.02)
    min_residual = min(cold.residuals.min(), warm.residuals.min())
    max_residual = max(cold.residuals.max(), warm.residuals.max())
    ax.set_ylim(min_residual * 0.9, max_residual * 1.1)

    cold_line, = ax.plot([], [], label="Cold start", color="#1f77b4")
    warm_line, = ax.plot([], [], label="Warm start", color="#ff7f0e")
    ax.legend()

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
                              output_path: Path) -> None:
    if not snapshots:
        raise ValueError("No snapshots available for animation")

    sorted_snaps = sorted(snapshots, key=lambda snap: snap.iteration)
    global_min = min(np.min(snap.field) for snap in sorted_snaps)
    global_max = max(np.max(snap.field) for snap in sorted_snaps)

    fig, (ax_field, ax_curve) = plt.subplots(1, 2, figsize=(10, 4))

    initial = sorted_snaps[0]
    im = ax_field.imshow(initial.field,
                         extent=[initial.x.min(), initial.x.max(), initial.y.min(), initial.y.max()],
                         origin="lower",
                         cmap="viridis",
                         vmin=global_min,
                         vmax=global_max,
                         aspect="auto")
    ax_field.set_title(f"Iteration {initial.iteration}")
    ax_field.set_xlabel("x (m)")
    ax_field.set_ylabel("y (m)")
    fig.colorbar(im, ax=ax_field, shrink=0.85, label="Az")

    ax_curve.set_title("Residual history")
    ax_curve.set_xlabel("Iteration")
    ax_curve.set_ylabel("Relative residual")
    ax_curve.set_yscale("log")
    ax_curve.set_xlim(0, history.iterations.max() * 1.02)
    ax_curve.set_ylim(history.residuals.min() * 0.9, history.residuals.max() * 1.1)

    progress_line, = ax_curve.plot([], [], color="#1f77b4")
    marker, = ax_curve.plot([], [], marker="o", color="#d62728")

    def update(frame: int):
        snapshot = sorted_snaps[frame]
        im.set_data(snapshot.field)
        im.set_extent([snapshot.x.min(), snapshot.x.max(), snapshot.y.min(), snapshot.y.max()])
        ax_field.set_title(f"Iteration {snapshot.iteration}")

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


def discover_frame_history(base_dir: Path, frame_index: int, warm: bool) -> Path:
    subdir = "warm" if warm else "cold"
    path = base_dir / "warm_start" / subdir / f"frame_{frame_index}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing history for frame {frame_index} ({'warm' if warm else 'cold'}): {path}")
    return path


def discover_snapshot_files(base_dir: Path) -> Tuple[HistorySeries, List[SnapshotRecord]]:
    history_path = base_dir / "progress_snapshots" / "history.csv"
    if not history_path.exists():
        raise FileNotFoundError(f"Snapshot history not found: {history_path}")
    history = load_history_csv(history_path)

    snapshot_dir = base_dir / "progress_snapshots" / "snapshots"
    snapshot_paths = sorted(snapshot_dir.glob("snapshot_iter_*.csv"))
    snapshots = [load_snapshot_csv(path) for path in snapshot_paths]
    return history, snapshots


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
    cold_path = discover_frame_history(dataset_root, args.frame, warm=False)
    warm_path = discover_frame_history(dataset_root, args.frame, warm=True)

    cold_history = load_history_csv(cold_path)
    warm_history = load_history_csv(warm_path)

    animations_dir = args.output_dir
    warm_animation_path = animations_dir / f"warm_start_frame_{args.frame}.gif"
    warm_image_path = animations_dir / f"warm_start_frame_{args.frame}.png"
    print(f"[viz] Creating warm-start animation at {warm_animation_path}")
    create_warm_start_animation(cold_history, warm_history, args.frame, warm_animation_path, warm_image_path)

    history, snapshots = discover_snapshot_files(dataset_root)
    snapshot_animation_path = animations_dir / "progress_snapshots.gif"
    print(f"[viz] Creating progress snapshot animation at {snapshot_animation_path}")
    create_snapshot_animation(history, snapshots, snapshot_animation_path)

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

